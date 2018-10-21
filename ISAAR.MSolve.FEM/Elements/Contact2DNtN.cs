using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Contact2DNtN : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[2] { DOFType.X, DOFType.Y };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double youngModulus;
        private IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public Contact2DNtN(double youngModulus)
        {
            this.youngModulus = youngModulus;
        }

        public Contact2DNtN(double youngModulus, IElementDOFEnumerator dofEnumerator)
            : this(youngModulus)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public int ID
        {
            get { return 1; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.TwoD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        private Vector CalculateNormalUnitVector(Element element)
        {
            double X1 = element.INodes[0].X;
            double Y1 = element.INodes[0].Y;
            double X2 = element.INodes[1].X;
            double Y2 = element.INodes[1].Y;
            Vector normalVector = new Vector(new double[] { X2 - X1, Y2 - Y1 });
            double normalVectorLength = normalVector.Norm;
            Vector normalUnitVec = new Vector(new double[] { normalVector[0] / normalVectorLength, normalVector[1] / normalVectorLength });
            return normalUnitVec;
        }

        private Matrix2D CalculatePositionMatrix()
        {
            Matrix2D aMatrix = new Matrix2D(new double[,]
                {
                    { -1,0,1,0},
                    {0,-1,0,1 }
                });
            return aMatrix;
        }

        private double CalculateNormalGap(Element element)
        {
            Matrix2D A = CalculatePositionMatrix();
            Matrix2D AT = A.Transpose();
            Vector n = CalculateNormalUnitVector(element);
            Vector AT_n = AT * n;
            double[] xupd = new double[] {
                Nodes[1].XCoordinate + DisplacementVector[0],
                Nodes[1].YCoordinate + DisplacementVector[1],
                Nodes[2].XCoordinate + DisplacementVector[2],
                Nodes[2].YCoordinate + DisplacementVector[3]
            };
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            double penetration = CalculateNormalGap();
            if (penetration <= 0)
            {
                double[] n = CalculateNormalUnitVector();
                double[,] A = CalculatePositionMatrix();
                double[,] AT = MatrixOperations.Transpose(A);
                double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
                double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
                double[,] AT_nxn_A = MatrixOperations.MatrixProduct(AT, nxn_A);
                double[,] globalStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
                return globalStiffnessMatrix;
            }
            else
            {
                double[,] globalStifnessMatrix = new double[4, 4];
                return globalStifnessMatrix;
            }
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);

            double totalMass = Density * SectionArea * L;

            return new Matrix2D(new double[,]
            {
                { totalMass/2,0,0,0},
                {0,totalMass/2,0,0 },
                {0,0,totalMass/2,0 },
                {0,0,0,totalMass/2 }
            });
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        #region IFiniteElement Members


        public bool MaterialModified
        {
            get { return false; }
        }

        public void ResetMaterialModified()
        {
        }

        #endregion

        #region IFiniteElement Members

        public void ClearMaterialState()
        {
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        #endregion

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] local_Displacements, double[] local_d_Displacements)
        {
            throw new Exception("Not implemented");
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new Exception("Not implemented");
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new Exception("Not implemented");
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }
    }
}
