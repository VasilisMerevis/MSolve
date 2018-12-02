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
        private double[] DisplacementVector { get; set; }
        double PenaltyFactor { get; set; }

        public Contact2DNtN(double youngModulus)
        {
            this.youngModulus = youngModulus;
            DisplacementVector = new double[4];
            PenaltyFactor = 100.0 * youngModulus;
        }

        public Contact2DNtN(double youngModulus, IElementDOFEnumerator dofEnumerator)
            : this(youngModulus)
        {
            this.dofEnumerator = dofEnumerator;
            DisplacementVector = new double[4];
            PenaltyFactor = 100.0 * youngModulus;
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

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
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
            Vector xupd = new Vector(new double[] {
                element.INodes[0].X + DisplacementVector[0],
                element.INodes[0].Y + DisplacementVector[1],
                element.INodes[1].X + DisplacementVector[2],
                element.INodes[1].Y + DisplacementVector[3]
            });
            double normalGap = xupd * AT_n;
            return normalGap;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            double penetration = CalculateNormalGap(element as Element);
            if (penetration <= 0)
            {
                Vector n = CalculateNormalUnitVector(element as Element);
                Matrix2D A = CalculatePositionMatrix();
                Matrix2D AT = A.Transpose();
                Matrix2D nxn = n.OuterProduct(n);
                Matrix2D nxn_A = nxn * A;
                Matrix2D AT_nxn_A = AT * nxn_A;
                AT_nxn_A.Scale(PenaltyFactor);
                IMatrix2D globalStiffnessMatrix = AT_nxn_A;
                return globalStiffnessMatrix;
            }
            else
            {
                IMatrix2D globalStifnessMatrix = new Matrix2D(new double[4, 4]);
                return globalStifnessMatrix;
            }
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            return null;
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            return null;
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
            return null;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            DisplacementVector = localDisplacements;
            double penetration = CalculateNormalGap(element);
            if (penetration <= 0)
            {
                Matrix2D A = CalculatePositionMatrix();
                Matrix2D AT = A.Transpose();
                Vector n = CalculateNormalUnitVector(element);
                Vector AT_n = AT * n;
                double ksi = CalculateNormalGap(element);
                Vector ksi_AT_n = ksi * AT_n;
                Vector e_ksi_AT_n = PenaltyFactor * ksi_AT_n;
                return e_ksi_AT_n.Data;
            }
            else
            {
                double[] internalGlobalForcesVector = new double[4];
                return internalGlobalForcesVector;
            }
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new Exception("Not implemented");
        }

        public void SaveMaterialState()
        {
            
        }
    }
}
