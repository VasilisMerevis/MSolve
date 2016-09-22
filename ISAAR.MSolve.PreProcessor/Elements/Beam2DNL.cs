using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    public class Beam2DNL : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.RotZ };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly IFiniteElementMaterial material;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public double SectionArea { get; set; }
        public double MomentOfInertia { get; set; }

        private Vector<double> internalLocalForcesVector = new Vector<double>(3);
        private Vector<double> localDisplacementVector = new Vector<double>(3);
        private double cosInitial, sinInitial, lengthInitial, betaAngleInitial;
        private double cosCurrent, sinCurrent, lengthCurrent, betaAngleCurrent;
        private double[] node1XYinitial, node2XYinitial, node1XYcurrent, node2XYcurrent;
        private int iter = 0;



        public Beam2DNL(IFiniteElementMaterial material)
        {
            this.material = material;
        }

        public Beam2DNL(IFiniteElementMaterial material, IFiniteElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        #region IElementType Members

        public int ID
        {
            get { return 1; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.TwoD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }
        
        private void GetInitialGeometricData(Element element)
        {
            node1XYinitial = new double[] { element.Nodes[0].X, element.Nodes[0].Y };
            node2XYinitial = new double[] { element.Nodes[1].X, element.Nodes[1].Y };
            lengthInitial = Math.Sqrt(Math.Pow((node2XYinitial[0] - node1XYinitial[0]), 2) + Math.Pow((node2XYinitial[1] - node1XYinitial[1]), 2));
            betaAngleInitial = Math.Atan2(node2XYinitial[1] - node1XYinitial[1], node2XYinitial[0] - node1XYinitial[0]);
            cosInitial = (node2XYinitial[0] - node1XYinitial[0]) / lengthInitial;
            sinInitial = (node2XYinitial[1] - node1XYinitial[1]) / lengthInitial;

            lengthCurrent = lengthInitial;
            sinCurrent = sinInitial;
            cosCurrent = cosInitial;
            betaAngleCurrent = betaAngleInitial;
        }

        private void GetCurrentGeometricalData(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            node1XYcurrent = new double[] {
                node1XYinitial[0]+localDisplacements[0]+localdDisplacements[0],
                node1XYinitial[1]+localDisplacements[1]+localdDisplacements[1]
                };
            node2XYcurrent = new double[] {
                node2XYinitial[0]+localDisplacements[3]+localdDisplacements[3],
                node2XYinitial[1]+localDisplacements[4]+localdDisplacements[4]
                };
            betaAngleCurrent = Math.Atan2(node2XYcurrent[1] - node1XYcurrent[1], node2XYcurrent[0] - node1XYcurrent[0]);

            lengthCurrent = Math.Sqrt(Math.Pow((node2XYcurrent[0] - node1XYcurrent[0]), 2) + Math.Pow((node2XYcurrent[1] - node1XYcurrent[1]), 2));
            cosCurrent = (node2XYcurrent[0] - node1XYcurrent[0]) / lengthCurrent;
            sinCurrent = (node2XYcurrent[1] - node1XYcurrent[1]) / lengthCurrent;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            //throw new Exception("Breakpoint here K");
            iter++;
            if (iter < 2)
            {
                GetInitialGeometricData(element);
            }
            else
            {
                //Console.WriteLine(iter);
            }
            
            //GetInitialGeometricData(element);

            double N = internalLocalForcesVector[0];
            double M1 = internalLocalForcesVector[1];
            double M2 = internalLocalForcesVector[2];
            double sinb = sinCurrent;
            double cosb = cosCurrent;
            double cosb2 = cosCurrent * cosCurrent;
            double sinb2 = sinCurrent * sinCurrent;
            double cosbsinb = cosCurrent * sinCurrent;
            double Lc = lengthCurrent;
            double L0 = lengthInitial;
            double Lc2 = lengthCurrent * lengthCurrent;
            double E = (material as ElasticMaterial).YoungModulus;
            double I = MomentOfInertia;
            double A = SectionArea;

            double[,] globalStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix[0, 0] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            globalStiffnessMatrix[0, 1] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            globalStiffnessMatrix[0, 2] = -6 * E * I * sinb / (L0 * Lc);
            globalStiffnessMatrix[0, 3] = -N * sinb2 / Lc - 12 * E * I * sinb2 / (L0 * Lc2) + 2 * (M2 + M1) * cosbsinb / Lc2 - A * E * cosb2 / L0;
            globalStiffnessMatrix[0, 4] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            globalStiffnessMatrix[0, 5] = -6 * E * I * sinb / (L0 * Lc);

            globalStiffnessMatrix[1, 0] = globalStiffnessMatrix[0, 1];
            globalStiffnessMatrix[1, 1] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[1, 2] = 6 * E * I * cosb / (L0 * Lc);
            globalStiffnessMatrix[1, 3] = (M2 + M1) * (sinb2 - cosb2) / Lc2 + N * cosbsinb / Lc + 12 * E * I * cosbsinb / (L0 * Lc2) - A * E * cosbsinb / L0;
            globalStiffnessMatrix[1, 4] = -A * E * sinb2 / L0 - 2 * (M2 + M1) * cosbsinb / Lc2 - N * cosb2 / Lc - 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[1, 5] = 6 * E * I * cosb / (L0 * Lc);

            globalStiffnessMatrix[2, 0] = globalStiffnessMatrix[0, 2];
            globalStiffnessMatrix[2, 1] = globalStiffnessMatrix[1, 2];
            globalStiffnessMatrix[2, 2] = 4 * E * I / L0;
            globalStiffnessMatrix[2, 3] = 6 * E * I * sinb / (L0 * Lc);
            globalStiffnessMatrix[2, 4] = -6 * E * I * cosb / (L0 * Lc);
            globalStiffnessMatrix[2, 5] = 2 * E * I / L0;

            globalStiffnessMatrix[3, 0] = globalStiffnessMatrix[0, 3];
            globalStiffnessMatrix[3, 1] = globalStiffnessMatrix[1, 3];
            globalStiffnessMatrix[3, 2] = globalStiffnessMatrix[2, 3];
            globalStiffnessMatrix[3, 3] = N * sinb2 / Lc + 12 * E * I * sinb2 / (L0 * Lc2) - 2 * (M2 + M1) * cosbsinb / Lc2 + A * E * cosb2 / L0;
            globalStiffnessMatrix[3, 4] = (M2 + M1) * (cosb2 - sinb2) / Lc2 - N * cosbsinb / Lc - 12 * E * I * cosbsinb / (L0 * Lc2) + A * E * cosbsinb / L0;
            globalStiffnessMatrix[3, 5] = 6 * E * I * sinb / (L0 * Lc);

            globalStiffnessMatrix[4, 0] = globalStiffnessMatrix[0, 4];
            globalStiffnessMatrix[4, 1] = globalStiffnessMatrix[1, 4];
            globalStiffnessMatrix[4, 2] = globalStiffnessMatrix[2, 4];
            globalStiffnessMatrix[4, 3] = globalStiffnessMatrix[3, 4];
            globalStiffnessMatrix[4, 4] = A * E * sinb2 / L0 + 2 * (M2 + M1) * cosbsinb / Lc2 + N * cosb2 / Lc + 12 * E * I * cosb2 / (L0 * Lc2);
            globalStiffnessMatrix[4, 5] = -6 * E * I * cosb / (L0 * Lc);

            globalStiffnessMatrix[5, 0] = globalStiffnessMatrix[0, 5];
            globalStiffnessMatrix[5, 1] = globalStiffnessMatrix[1, 5];
            globalStiffnessMatrix[5, 2] = globalStiffnessMatrix[2, 5];
            globalStiffnessMatrix[5, 3] = globalStiffnessMatrix[3, 5];
            globalStiffnessMatrix[5, 4] = globalStiffnessMatrix[4, 5];
            globalStiffnessMatrix[5, 5] = 4 * E * I / L0;

            IMatrix2D<double> iGlobalStiffnessMatrix = new Matrix2D<double>(globalStiffnessMatrix);

            return dofEnumerator.GetTransformedMatrix(iGlobalStiffnessMatrix);
        }

        public IMatrix2D<double> MassMatrix(Element element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double L2 = L * L;
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double c2 = c * c;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
            double s2 = s * s;
            double dAL420 = Density * SectionArea * L / 420;

            double totalMass = Density * SectionArea * L;
            double totalMassOfDiagonalTerms = 2 * dAL420 * (140 * c2 + 156 * s2) + 2 * dAL420 * (140 * s2 + 156 * c2);
            double scale = totalMass / totalMassOfDiagonalTerms;

            return new SymmetricMatrix2D<double>(new double[] { dAL420*(140*c2+156*s2)*scale, 0, 0, 0, 0, 0,
                dAL420*(140*s2+156*c2)*scale, 0, 0, 0, 0,
                0, 0, 0, 0,
                dAL420*(140*c2+156*s2)*scale, 0, 0,
                dAL420*(140*s2+156*c2)*scale, 0,
                0 });
        }

        public IMatrix2D<double> DampingMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return null;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        //public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        //{
        //    throw new NotImplementedException();
        //}
        public void CalculateLocalDisplacementVector(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            this.localDisplacementVector[0] = lengthCurrent - lengthInitial;
            this.localDisplacementVector[1] = (localDisplacements[2] + localdDisplacements[2]) - betaAngleCurrent + betaAngleInitial;
            this.localDisplacementVector[2] = (localDisplacements[5] + localdDisplacements[5]) - betaAngleCurrent + betaAngleInitial;
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            //throw new Exception("Breakpoint here F");
            double E = (material as ElasticMaterial).YoungModulus;
            double I = MomentOfInertia;
            double A = SectionArea;

            GetCurrentGeometricalData(element, localDisplacements, localdDisplacements);

            
            Matrix2D<double> Dmatrix = new Matrix2D<double>(new[,]
            {
                { E * A / lengthInitial, 0, 0 },
                { 0 , 4 * E * I / lengthInitial, 2 * E * I / lengthInitial },
                { 0 , 2 * E * I / lengthInitial, 4 * E * I / lengthInitial }
            });

            Matrix2D<double> Bmatrix = new Matrix2D<double>(new[,]
            {
                { -cosCurrent, -sinCurrent, 0, cosCurrent, sinCurrent, 0 },
                { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 1, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 0 },
                { -sinCurrent / lengthCurrent, cosCurrent / lengthCurrent, 0, sinCurrent / lengthCurrent, -cosCurrent / lengthCurrent, 1 }
            });

            CalculateLocalDisplacementVector(element, localDisplacements, localdDisplacements);
            Vector<double> internalLocalForcesVector = Dmatrix * localDisplacementVector;
            Vector<double> internalGlobalForcesVector = Bmatrix.Transpose() * internalLocalForcesVector;
            double[] internalGlobalForces = internalGlobalForcesVector.Data;// as double[];
           
            return internalGlobalForces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector<double> accelerations = new Vector<double>(6);
            IMatrix2D<double> massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            double[] forces = new double[6];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void SaveMaterialState()
        {
            //do nothing
        }

        #endregion

        #region IStructuralFiniteElement Members

        public IFiniteElementMaterial Material
        {
            get { return material; }
        }

        #endregion

        #region IFiniteElement Members


        public bool MaterialModified
        {
            get { return material.Modified; }
        }

        public void ResetMaterialModified()
        {
            material.ResetModified();
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
    }
}
