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
    public class Contact3DNtS : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        private readonly IFiniteElementMaterial material;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }


        private bool isInitializedK = false;
        private bool isInitializedF = false;

        private IMatrix2D<double> A;
        private Dictionary<int, IMatrix2D<double>> dA;
        private Dictionary<int, IMatrix2D<double>> ddA;
        private Vector<double> normalVector;
        private double penaltyFactor;
        private Vector<double> xInitial, xUpdated;
        private bool contactStatus;
        private Dictionary<int, double> N, dN, ddN;
        private Dictionary<int, Vector<double>> dRho, ddRho;
        private Matrix2D<double> metricTensor, inverseMetricTensor;
        private double detm, ksi3Penetration;

        public Contact3DNtS(IFiniteElementMaterial3D material)
        {
            this.material = material;
            this.penaltyFactor = material.YoungModulus;
        }

        public Contact3DNtS(IFiniteElementMaterial3D material, IFiniteElementDOFEnumerator dofEnumerator)
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
            normalVector = (1/Math.Sqrt(detm)) * Vector<double>.CrossProductInR3(dRho[1], dRho[2]);
        }

        private Vector<double> CalclulateDeltaKsi (Vector<double> xUpdated, Dictionary<int, Vector<double>> dRho, Dictionary<int, Vector<double>> ddRho, Matrix2D<double> A, Matrix2D<double> metricTensor, double detm)
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
                Console.WriteLine(deltaKsi.Norm);
                if (deltaKsi.Norm<0.01 && i<100)
                {
                    return ksiVector;
                }
                else if (i==100)
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

        private bool CheckContactStatus(Vector<double> xUpdated)
        {
            Vector<double> ksi = CPP();
            if (Math.Abs(ksi[0])>=1.05 && Math.Abs(ksi[1])>=1.05 )
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
            if (ksi3 >= 0.0)
            {
                Console.WriteLine("No penetration. No contact occurs");
                contactStatus = false;
                return contactStatus;
            }
            contactStatus = true;
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
            Matrix2D<double> mainPart = (posA.Transpose()*(normalVector ^ normalVector)*posA);
            mainPart.Scale(penaltyFactor);
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
            Console.WriteLine("InternalForcesCalculation");
            GetCurrentPosition(element, localDisplacements);
            if (ksi3Penetration>0.0)
            {
                double[] internalForcesVector = new double[15];
                return internalForcesVector;
            }
            
            Matrix2D<double> posA = (Matrix2D<double>)A;
            Vector<double> AT_n = posA.Transpose() * normalVector;
            double en_ksi3 = -penaltyFactor * ksi3Penetration;
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
