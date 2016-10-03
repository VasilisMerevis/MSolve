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

        private double cosInitial, sinInitial, lengthInitial, betaAngleInitial;
        private double cosCurrent, sinCurrent, lengthCurrent, betaAngleCurrent;
        private double[] node1XYInitial, node2XYInitial, node1XYCurrent, node2XYCurrent;
        private double[] node1GlobalDisplacementVector, node2GlobalDisplacementVector;
        private int iter = 0;

        private int isInitialized = 0;



        public Beam2DNL(IFiniteElementMaterial material)
        {
            this.material = material;
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];

        }

        public Beam2DNL(IFiniteElementMaterial material, IFiniteElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
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
            node1XYInitial = new double[] { element.Nodes[0].X, element.Nodes[0].Y };
            node2XYInitial = new double[] { element.Nodes[1].X, element.Nodes[1].Y };
            lengthInitial = Math.Sqrt(Math.Pow(node2XYInitial[0] - node1XYInitial[0], 2) + Math.Pow(node2XYInitial[1] - node1XYInitial[1], 2));
            betaAngleInitial = Math.Atan2(node2XYInitial[1] - node1XYInitial[1], node2XYInitial[0] - node1XYInitial[0]);
            cosInitial = (node2XYInitial[0] - node1XYInitial[0]) / lengthInitial;
            sinInitial = (node2XYInitial[1] - node1XYInitial[1]) / lengthInitial;

            lengthCurrent = lengthInitial;
            sinCurrent = sinInitial;
            cosCurrent = cosInitial;
            betaAngleCurrent = betaAngleInitial;
        }

        private void GetCurrentGeometricalData()
        {
            node1XYCurrent = new[]
            {
                node1XYInitial[0] + node1GlobalDisplacementVector[0],
                node1XYInitial[1] + node1GlobalDisplacementVector[1]
            };

            node2XYCurrent = new[]
            {
                node2XYInitial[0] + node2GlobalDisplacementVector[0],
                node2XYInitial[1] + node2GlobalDisplacementVector[1]
            };
            betaAngleCurrent = Math.Atan2(node2XYCurrent[1] - node1XYCurrent[1], node2XYCurrent[0] - node1XYCurrent[0]);
            
            lengthCurrent = Math.Sqrt(Math.Pow(node2XYCurrent[0] - node1XYCurrent[0], 2) + Math.Pow(node2XYCurrent[1] - node1XYCurrent[1], 2));
            cosCurrent = (node2XYCurrent[0] - node1XYCurrent[0]) / lengthCurrent;
            sinCurrent = (node2XYCurrent[1] - node1XYCurrent[1]) / lengthCurrent;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            //throw new Exception("K");
            if (this.isInitialized == 0)
            {
                this.GetInitialGeometricData(element);
                this.isInitialized = 1;
            }
            else
            {
                Console.WriteLine("else");
            }

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
            double N = A * E * (lengthCurrent - lengthInitial) / lengthInitial;
            double M1 = 2 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 4 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double M2 = 4 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 2 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;

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

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localDeltaDisplacements)
        {
            //throw new Exception("F");

            node1GlobalDisplacementVector = new double[] { localDisplacements[0] + localDeltaDisplacements[0], localDisplacements[1] + localDeltaDisplacements[1], localDisplacements[2] + localDeltaDisplacements[2] };
            node2GlobalDisplacementVector = new double[] { localDisplacements[3] + localDeltaDisplacements[3], localDisplacements[4] + localDeltaDisplacements[4], localDisplacements[5] + localDeltaDisplacements[5] };
            double[] totalDisplacementVector = new double[localDisplacements.Length];
            for(int i = 0; i < localDisplacements.Length; ++i)
            {
                totalDisplacementVector[i] = localDisplacements[i] + localDeltaDisplacements[i];
            }
            
            GetCurrentGeometricalData();
            double E = (material as ElasticMaterial).YoungModulus;
            double I = MomentOfInertia;
            double A = SectionArea;
            double N = A * E * (lengthCurrent - lengthInitial) / lengthInitial;
            double M1 = 2 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 4 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double M2 = 4 * E * I * (node2GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial + 2 * E * I * (node1GlobalDisplacementVector[2] - betaAngleCurrent + betaAngleInitial) / lengthInitial;
            double[] intforceV = new double[6];
            intforceV[0] = -(M2 * sinCurrent / lengthCurrent) - (M1 * sinCurrent / lengthCurrent) - N * cosCurrent;
            intforceV[1] = (M2 * cosCurrent / lengthCurrent) + (M1 * cosCurrent / lengthCurrent) - N * sinCurrent;
            intforceV[2] = M1;
            intforceV[3] = (M2 * sinCurrent / lengthCurrent) + (M1 * sinCurrent / lengthCurrent) + N * cosCurrent;
            intforceV[4] = -(M2 * cosCurrent / lengthCurrent) - (M1 * cosCurrent / lengthCurrent) + N * sinCurrent;
            intforceV[5] = M2;

            return intforceV;
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
