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
    public class Contact3DNtSFr : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        private readonly IFiniteElementMaterial material;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public bool RotationalStiffness { get; set; }


        private bool isInitializedK = false;
        private bool isInitializedF = false;

        private IMatrix2D<double> A;
        private Dictionary<int, IMatrix2D<double>> dA;
        private Dictionary<int, IMatrix2D<double>> ddA;
        private Vector<double> normalVector;
        private double penaltyFactor, tangentPenaltyFactor;
        private Vector<double> xInitial, xUpdated;
        private bool contactStatus;
        private Dictionary<int, double> N, dN, ddN;
        private Dictionary<int, Vector<double>> dRho, ddRho, dRhoPrevious;
        private Matrix2D<double> metricTensor, inverseMetricTensor;
        private double detm, ksi3Penetration, frictionCoef;
        private Vector<double> Tr, TrPrevious;
        private Vector<double> rho, rhoPrevious;

        public Contact3DNtSFr(IFiniteElementMaterial3D material)
        {
            this.material = material;
            this.penaltyFactor = material.YoungModulus * 100000.0;
            RotationalStiffness = true;
        }

        public Contact3DNtSFr(IFiniteElementMaterial3D material, IFiniteElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public int ID
        {
            get { return 1; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        private void GetInitialPosition(Element element)
        {
            double[] node1XYZInitial = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
            double[] node2XYZInitial = new double[] { element.Nodes[1].X, element.Nodes[1].Y, element.Nodes[1].Z };
            double[] node3XYZInitial = new double[] { element.Nodes[2].X, element.Nodes[2].Y, element.Nodes[2].Z };
            double[] node4XYZInitial = new double[] { element.Nodes[3].X, element.Nodes[3].Y, element.Nodes[3].Z };
            double[] node5XYZInitial = new double[] { element.Nodes[4].X, element.Nodes[4].Y, element.Nodes[4].Z };
            double[] xVector = new double[15];
            node1XYZInitial.CopyTo(xVector, 0);
            node2XYZInitial.CopyTo(xVector, node1XYZInitial.Length);
            node3XYZInitial.CopyTo(xVector, node1XYZInitial.Length + node2XYZInitial.Length);
            node4XYZInitial.CopyTo(xVector, node1XYZInitial.Length + node2XYZInitial.Length + node3XYZInitial.Length);
            node5XYZInitial.CopyTo(xVector, node1XYZInitial.Length + node2XYZInitial.Length + node3XYZInitial.Length + node4XYZInitial.Length);
            xInitial = new Vector<double>(xVector);
            xUpdated = new Vector<double>(xVector);
        }

        private void GetCurrentPosition(Element element, double[] displacements)
        {
            Vector<double> displacementVector = new Vector<double>(displacements);
            double[] xUpdatedVetor = xInitial + displacementVector;
            xUpdated = new Vector<double>(xUpdatedVetor);
        }

        private void CalculateShapeFunctions(double ksi1, double ksi2)
        {
            N = new Dictionary<int, double>();
            dN = new Dictionary<int, double>();
            ddN = new Dictionary<int, double>();
            //Shape functions
            N[1] = (1.0 / 4.0) * (1.0 - ksi1) * (1.0 - ksi2);
            N[2] = (1.0 / 4.0) * (1.0 + ksi1) * (1.0 - ksi2);
            N[3] = (1.0 / 4.0) * (1.0 + ksi1) * (1.0 + ksi2);
            N[4] = (1.0 / 4.0) * (1.0 - ksi1) * (1.0 + ksi2);

            //First order derivatives of shape functions
            dN[11] = -(1.0 / 4.0) * (1.0 - ksi2);
            dN[12] = -(1.0 / 4.0) * (1.0 - ksi1);
            dN[21] = (1.0 / 4.0) * (1.0 - ksi2);
            dN[22] = -(1.0 / 4.0) * (1.0 + ksi1);
            dN[31] = (1.0 / 4.0) * (1.0 + ksi2);
            dN[32] = (1.0 / 4.0) * (1.0 + ksi1);
            dN[41] = -(1.0 / 4.0) * (1.0 + ksi2);
            dN[42] = (1.0 / 4.0) * (1.0 - ksi1);

            //Second order derivatives of shape functions
            ddN[112] = 1.0 / 4.0;
            ddN[212] = -1.0 / 4.0;
            ddN[312] = 1.0 / 4.0;
            ddN[412] = -1.0 / 4.0;
        }

        /// <summary>
        /// Returns position matrix A as well as its first and secoond order derivatives dA and ddA respectively
        /// </summary>
        /// <param name="N">Shape functions</param>
        /// <param name="dN">Shape functions first order derivatives</param>
        /// <param name="ddN">Shape functions second order derivatives</param>
        /// <returns>Returns position matrix A as well as its first and secoond order derivatives dA and ddA respectively</returns>
        private void PositionMatrices()
        {
            dA = new Dictionary<int, IMatrix2D<double>>();
            ddA = new Dictionary<int, IMatrix2D<double>>();
            //Position matrix A
            double[,] AMatrix = new double[,]
            {
                {-N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1.0, 0, 0},
                {0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1.0, 0 },
                {0, 0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1.0 }
            };
            A = new Matrix2D<double>(AMatrix);

            //First order derivatives of position matrix
            double[,] dA1Matrix = new double[,]
            {
                {-dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0, 0, 0},
                {0, -dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0, 0 },
                {0, 0, -dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0 }
            };

            double[,] dA2Matrix = new double[,]
            {
                {-dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0, 0, 0},
                {0, -dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0, 0 },
                {0, 0, -dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0 }
            };

            dA[1] = new Matrix2D<double>(dA1Matrix);
            dA[2] = new Matrix2D<double>(dA2Matrix);

            //Second order derivatives of position matrix
            double[,] ddA12 = new double[,]
            {
                {-ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0, 0},
                {0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0 },
                {0, 0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0 }
            };

            ddA[12] = new Matrix2D<double>(ddA12);
        }

        private void SurfaceGeometry()
        {
            dRho = new Dictionary<int, Vector<double>>();
            ddRho = new Dictionary<int, Vector<double>>();
            //Surface vectors
            dRho[1] = -1.0 * ((Matrix2D<double>)dA[1] * (Vector<double>)xUpdated);
            dRho[2] = -1.0 * ((Matrix2D<double>)dA[2] * (Vector<double>)xUpdated);

            //Surface vectors derivatives
            ddRho[12] = -1.0 * ((Matrix2D<double>)ddA[12] * (Vector<double>)xUpdated);

            //Metric properties
            metricTensor = new Matrix2D<double>(new double[,]
                { { dRho[1]*dRho[1], dRho[1]*dRho[2]}, {dRho[2]*dRho[1], dRho[2]*dRho[2] } }
                );

            detm = metricTensor[0, 0] * metricTensor[1, 1] - metricTensor[1, 0] * metricTensor[0, 1];

            inverseMetricTensor = new Matrix2D<double>(new double[,] {
                { metricTensor[1,1]/detm, -metricTensor[0,1]/detm},
                { -metricTensor[1,0]/detm, metricTensor[0,0]/detm }
            });

            //Normal vector
            normalVector = (1 / Math.Sqrt(detm)) * Vector<double>.CrossProductInR3(dRho[1], dRho[2]);
        }

        private Vector<double> CalclulateDeltaKsi(Vector<double> xUpdated, Dictionary<int, Vector<double>> dRho, Dictionary<int, Vector<double>> ddRho, Matrix2D<double> A, Matrix2D<double> metricTensor, double detm)
        {
            double f1 = dRho[1].DotProduct(A * xUpdated);
            double f2 = dRho[2].DotProduct(A * xUpdated);
            double e = ddRho[12].DotProduct(A * xUpdated);
            Vector<double> f = new Vector<double>(new double[] { f1, f2 });
            double[,] matrixM = new double[,] { { metricTensor[1, 1], e - metricTensor[0, 1] }, { e - metricTensor[1, 0], metricTensor[0, 0] } };
            Matrix2D<double> interMatrix = new Matrix2D<double>(matrixM);

            double coef = 1.0 / (detm - Math.Pow(e, 2) + 2.0 * e * metricTensor[0, 1]);
            Vector<double> deltaKsi = coef * (interMatrix * f);
            return deltaKsi;
        }

        private Vector<double> CPP()
        {
            double ksi1 = 0.0;
            double ksi2 = 0.0;
            Vector<double> ksiVector = new Vector<double>(new double[] { ksi1, ksi2 });
            for (int i = 0; i <= 100; i++)
            {
                CalculateShapeFunctions(ksiVector.Data[0], ksiVector.Data[1]);
                PositionMatrices();
                SurfaceGeometry();
                Vector<double> deltaKsi = CalclulateDeltaKsi(xUpdated, dRho, ddRho, (Matrix2D<double>)A, metricTensor, detm);
                //Console.WriteLine(deltaKsi.Norm);
                if (deltaKsi.Norm < 0.01 && i < 100)
                {
                    return ksiVector;
                }
                else if (i == 100)
                {

                    throw new Exception("Contact did not converge. ");
                    break;
                }
                else
                {
                    double[] ksi = ksiVector + deltaKsi;
                    ksiVector = new Vector<double>(ksi);
                }

            }
            return ksiVector;
        }

        private double Penetration(Matrix2D<double> A, Vector<double> xUpdated, Vector<double> n)
        {
            double ksi3 = xUpdated.DotProduct(A.Transpose() * n);
            return ksi3;
        }

        private Vector<double> CalculateCurrentTangentialTraction(double[] displacementVector)
        {
            double[] tangentialTraction = new double[2];
            double b1 = dRhoPrevious[1] * dRho[1];
            double b2 = dRhoPrevious[2] * dRho[1];
            double c1 = dRhoPrevious[1] * dRho[2];
            double c2 = dRhoPrevious[2] * dRho[2];
            double[,] m = metricTensor.Data;
           
            double[,] ArMatrix = new double[,]
            {
                {-N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0},
                {0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0},
                {0, 0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4]}
            };
            Matrix2D<double> Ar = new Matrix2D<double>(ArMatrix);
            double[] masterDisplacementVector = new double[12];
            double[] projectionPointCoordinatesVector = new double[12];
            for (int i = 0; i < masterDisplacementVector.Length; i++)
            {
                masterDisplacementVector[i] = displacementVector[i];
                projectionPointCoordinatesVector[i] = xUpdated[i];
            }
            Vector<double> u = new Vector<double>(masterDisplacementVector);
            Vector<double> xcUpdated = new Vector<double>(projectionPointCoordinatesVector);
            rhoPrevious = rho;
            rho = Ar * xcUpdated;
            
            Vector<double> rhoPrevius_plus_u = new Vector<double>(rhoPrevious + u);
            Vector<double> deltaRho = new Vector<double>(rho - rhoPrevius_plus_u);

            double d1 = deltaRho * dRho[1];
            double d2 = deltaRho * dRho[2];
            double et = tangentPenaltyFactor;
            tangentialTraction[0] = (TrPrevious[0] * m[0, 0] + TrPrevious[1] * m[1, 0]) * b1 + (TrPrevious[0] * m[0, 1] + TrPrevious[1] * m[1, 1]) * b2 - (et * d1);
            tangentialTraction[1] = (TrPrevious[0] * m[0, 0] + TrPrevious[1] * m[1, 0]) * c1 + (TrPrevious[0] * m[0, 1] + TrPrevious[1] * m[1, 1]) * c2 - (et * d2);
            Vector<double> Tr = new Vector<double>(tangentialTraction);
            return Tr;
        }

        private double CalculateCoulombFriction()
        {
            double phi;
            double normalForce = penaltyFactor * ksi3Penetration;
            double[,] m = metricTensor.Data;
            double sqRoot = Tr[0] * Tr[0] * m[0, 0] + Tr[0] * Tr[1] * m[1, 2] + Tr[1] * Tr[0] * m[1, 0] + Tr[1] * Tr[1] * m[1, 1];
            phi = Math.Sqrt(sqRoot) - (frictionCoef * normalForce);
            return phi;
        }

        private bool CheckContactStatus(Vector<double> xUpdated)
        {
            Vector<double> ksi = CPP();
            if (Math.Abs(ksi[0]) >= 1.05 && Math.Abs(ksi[1]) >= 1.05)
            {
                Console.WriteLine("Projection point out of surface. No contact occurs");
                contactStatus = false;
                return contactStatus;
            }
            CalculateShapeFunctions(ksi[0], ksi[1]);
            PositionMatrices();
            SurfaceGeometry();
            double ksi3 = Penetration((Matrix2D<double>)A, xUpdated, normalVector);
            ksi3Penetration = ksi3;
            if (ksi3 > 0.0)
            {
                Console.WriteLine("No penetration. No contact occurs");
                contactStatus = false;
                return contactStatus;
            }
            contactStatus = true;
            Console.WriteLine("Contact detected");
            return contactStatus;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            IMatrix2D<double> stiffnessMatrix;
            if (this.isInitializedK == false)
            {
                this.GetInitialPosition(element);
                this.isInitializedK = true;
            }

            bool activeContact = CheckContactStatus(xUpdated);
            if (activeContact == false)
            {
                stiffnessMatrix = new Matrix2D<double>(new double[15, 15]);
                return stiffnessMatrix;
            }

            Matrix2D<double> posA = (Matrix2D<double>)A;
            Matrix2D<double> mainPart = (posA.Transpose() * (normalVector ^ normalVector) * posA);
            mainPart.Scale(penaltyFactor);
            if (RotationalStiffness == true)
            {
                Matrix2D<double> posdA1 = (Matrix2D<double>)dA[1];
                Matrix2D<double> posdA2 = (Matrix2D<double>)dA[2];
                Matrix2D<double> part1 = (posdA1.Transpose() * (normalVector ^ dRho[1]) * posA) + (posA.Transpose() * (dRho[1] ^ normalVector) * posdA1);
                Matrix2D<double> part2 = (posdA1.Transpose() * (normalVector ^ dRho[2]) * posA) + (posA.Transpose() * (dRho[1] ^ normalVector) * posdA2);
                Matrix2D<double> part3 = (posdA2.Transpose() * (normalVector ^ dRho[1]) * posA) + (posA.Transpose() * (dRho[2] ^ normalVector) * posdA1);
                Matrix2D<double> part4 = (posdA2.Transpose() * (normalVector ^ dRho[2]) * posA) + (posA.Transpose() * (dRho[2] ^ normalVector) * posdA2);

                part1.Scale(penaltyFactor * ksi3Penetration * metricTensor[0, 0]);
                part2.Scale(penaltyFactor * ksi3Penetration * metricTensor[1, 0]);
                part3.Scale(penaltyFactor * ksi3Penetration * metricTensor[0, 1]);
                part4.Scale(penaltyFactor * ksi3Penetration * metricTensor[1, 1]);
                Matrix2D<double> rotationalPart = part1 + part2 + part3 + part4;
                stiffnessMatrix = mainPart + rotationalPart;
                return stiffnessMatrix;
            }


            //Sticking case
            Matrix2D<double> stickTangentialPart;
            Matrix2D<double> K2Sa;
            Matrix2D<double> K2Sb;
            Matrix2D<double> stickPart1 = new Matrix2D<double>(new double[15, 15]);
            Matrix2D<double> stickPart2 = new Matrix2D<double>(new double[15, 15]);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    Matrix2D<double> Apos = (Matrix2D<double>)A;
                    Matrix2D<double> K1S = Apos.Transpose() * (dRho[i] ^ dRho[j]) * Apos;
                    K1S.Scale(metricTensor[i, j]);
                    stickPart1 = stickPart1 + K1S; 
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            Matrix2D<double> dApos = (Matrix2D<double>)dA[j];
                            
                            K2Sa = (Apos.Transpose() * (dRho[k] ^ dRho[l]) * (Matrix2D<double>)dA[j]);
                            K2Sa.Scale(Tr[i] * metricTensor[i, l] * metricTensor[j, k]);
                            K2Sb = (dApos.Transpose() * (dRho[k] ^ dRho[l]) * (Matrix2D<double>)A);
                            K2Sb.Scale(Tr[i] * metricTensor[i, k] * metricTensor[j, l]);
                            stickPart2 = stickPart2 + K2Sa + K2Sb;
                        }
                    }
                }
            }
            stickTangentialPart = stickPart1 + stickPart2;
            //end of stickPart

            stiffnessMatrix = mainPart;
            return stiffnessMatrix;
        }

        public IMatrix2D<double> MassMatrix(Element element)
        {
            throw new NotImplementedException();
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
            //Console.WriteLine("InternalForcesCalculation");
            GetCurrentPosition(element, localDisplacements);
            bool activeContact = CheckContactStatus(xUpdated);
            Console.WriteLine("Displacement of master node 1 is: {0}", localDisplacements[1]);
            Console.WriteLine("Displacement of master node 2 is: {0}", localDisplacements[4]);
            Console.WriteLine("Displacement of master node 3 is: {0}", localDisplacements[7]);
            Console.WriteLine("Displacement of master node 4 is: {0}", localDisplacements[10]);
            Console.WriteLine("Displacement of slave node is: {0}", localDisplacements[13]);
            Console.WriteLine("Penetration is: {0}", ksi3Penetration);

            if (ksi3Penetration > 0.0)
            {
                double[] internalForcesVector = new double[15];
                return internalForcesVector;
            }
            //Console.WriteLine(ksi3Penetration);

            Matrix2D<double> posA = (Matrix2D<double>)A;
            Vector<double> AT_n = posA.Transpose() * normalVector;
            double en_ksi3 = penaltyFactor * ksi3Penetration;
            Vector<double> internalForce = en_ksi3 * AT_n;
            return internalForce.Data;
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

            double[] forces = new double[15];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void SaveMaterialState()
        {

        }

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
