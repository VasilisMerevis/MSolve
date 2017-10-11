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
    class Contact3DNtS : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        private readonly IFiniteElementMaterial material;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }

        private double[] node1XYInitial, node2XYInitial, node1XYCurrent, node2XYCurrent;
        private double[] node1GlobalDisplacementVector, node2GlobalDisplacementVector;

        private bool isInitializedK = false;
        private bool isInitializedF = false;

        private IMatrix2D<double> A;
        private Dictionary<int, IMatrix2D<double>> dA;
        private Dictionary<int, IMatrix2D<double>> ddA;
        private Vector<double> normalVector;
        private double penaltyFactor;

        public Contact3DNtS(IFiniteElementMaterial3D material)
        {
            this.material = material;
            this.node1GlobalDisplacementVector = new double[3];
            this.node2GlobalDisplacementVector = new double[3];
        }

        public Contact3DNtS(IFiniteElementMaterial3D material, IFiniteElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }

        private Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>> CalculateShapeFunctions(double ksi1, double ksi2)
        {
            double N1 = 1 / 4 * (1 - ksi1) * (1 - ksi2);
            double N2 = 1 / 4 * (1 + ksi1) * (1 - ksi2);
            double N3 = 1 / 4 * (1 + ksi1) * (1 + ksi2);
            double N4 = 1 / 4 * (1 - ksi1) * (1 + ksi2);

            Dictionary<int, double> N = new Dictionary<int, double>();
            N[1] = N1;
            N[2] = N2;
            N[3] = N3;
            N[4] = N4;

            double dN11 = -1 / 4 * (1 - ksi2);
            double dN12 = -1 / 4 * (1 - ksi1);
            double dN21 = 1 / 4 * (1 - ksi2);
            double dN22 = -1 / 4 * (1 + ksi1);
            double dN31 = 1 / 4 * (1 + ksi2);
            double dN32 = 1 / 4 * (1 + ksi1);
            double dN41 = -1 / 4 * (1 + ksi2);
            double dN42 = 1 / 4 * (1 - ksi1);

            Dictionary<int, double> dN = new Dictionary<int, double>();
            dN[11] = dN11;
            dN[12] = dN12;
            dN[21] = dN21;
            dN[22] = dN22;
            dN[31] = dN31;
            dN[32] = dN32;
            dN[41] = dN41;
            dN[42] = dN42;

            double ddN112 = 1 / 4;
            double ddN212 = -1 / 4;
            double ddN312 = 1 / 4;
            double ddN412 = -1 / 4;

            Dictionary<int, double> ddN = new Dictionary<int, double>();
            ddN[112] = ddN112;
            ddN[212] = ddN212;
            ddN[312] = ddN312;
            ddN[412] = ddN412;

            Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>> shapeFunctions = new Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>>(N, dN, ddN);
            return shapeFunctions;
        }

        /// <summary>
        /// Returns position matrix A as well as its first and secoond order derivatives dA and ddA respectively
        /// </summary>
        /// <param name="N">Shape functions</param>
        /// <param name="dN">Shape functions first order derivatives</param>
        /// <param name="ddN">Shape functions second order derivatives</param>
        /// <returns>Returns position matrix A as well as its first and secoond order derivatives dA and ddA respectively</returns>
        private Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>> PositionMatrices(Dictionary<int, double> N, Dictionary<int, double> dN, Dictionary<int, double> ddN)
        {
            double[,] A = new double[,]
            {
                {-N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1, 0, 0},
                {0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1, 0 },
                {0, 0, -N[1], 0, 0, -N[2], 0, 0, -N[3], 0, 0, -N[4], 0, 0, 1 }
            };
            IMatrix2D<double> AMatrix = new Matrix2D<double>(A);

            double[,] dA1 = new double[,]
            {
                {-dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0, 0, 0},
                {0, -dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0, 0 },
                {0, 0, -dN[11], 0, 0, -dN[21], 0, 0, -dN[31], 0, 0, -dN[41], 0, 0, 0 }
            };
            IMatrix2D<double> dA1Matrix = new Matrix2D<double>(dA1);

            double[,] dA2 = new double[,]
            {
                {-dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0, 0, 0},
                {0, -dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0, 0 },
                {0, 0, -dN[12], 0, 0, -dN[22], 0, 0, -dN[32], 0, 0, -dN[42], 0, 0, 0 }
            };
            IMatrix2D<double> dA2Matrix = new Matrix2D<double>(dA2);

            Dictionary<int, IMatrix2D<double>> dA = new Dictionary<int, IMatrix2D<double>>();
            dA[1] = dA1Matrix;
            dA[2] = dA2Matrix;

            double[,] ddA12 = new double[,]
            {
                {-ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0, 0},
                {0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0 },
                {0, 0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0 }
            };
            IMatrix2D<double> ddA12Matrix = new Matrix2D<double>(ddA12);

            Dictionary<int, IMatrix2D<double>> ddA = new Dictionary<int, IMatrix2D<double>>();
            ddA[12] = ddA12Matrix;

            Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>> positionMatrices = 
                new Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>>(AMatrix, dA, ddA);
            return positionMatrices;
        }

        private Tuple<Dictionary<int, Vector<double>>, Dictionary<int, Vector<double>>, Vector<double>, Matrix2D<double>, double> SurfaceGeometry(IVector<double> xUpdated, Dictionary<int, IMatrix2D<double>> dA, Dictionary<int, IMatrix2D<double>> ddA)
        {
            Vector<double> dRho1 = -1.0 * ((Matrix2D<double>)dA[1] * (Vector<double>)xUpdated);
            Vector<double> dRho2 = -1.0 * ((Matrix2D<double>)dA[2] * (Vector<double>)xUpdated);
            Dictionary<int, Vector<double>> dRho = new Dictionary<int, Vector<double>>();
            dRho[1] = dRho1;
            dRho[2] = dRho2;

            Vector<double> ddRho12 = -1.0 * ((Matrix2D<double>)ddA[12] * (Vector<double>)xUpdated);
            Dictionary<int, Vector<double>> ddRho = new Dictionary<int, Vector<double>>();
            ddRho[12] = ddRho12;

            Matrix2D<double> metricTensor = new Matrix2D<double>(new double[,]
                { { dRho1*dRho1, dRho1*dRho2}, {dRho2*dRho1, dRho2*dRho2 } }
                );

            double detm = metricTensor[0, 0] * metricTensor[1, 1] - metricTensor[1, 0] * metricTensor[0, 1];

            Matrix2D<double> inverseMetricTensor = new Matrix2D<double>(new double[,] {
                { metricTensor[1,1]/detm, -metricTensor[0,1]/detm},
                { -metricTensor[1,0]/detm, metricTensor[0,0]/detm }
            });

            Vector<double> normalVector = (1/Math.Sqrt(detm)) * Vector<double>.CrossProductInR3(dRho1, dRho2);

            return new Tuple<Dictionary<int, Vector<double>>, Dictionary<int, Vector<double>>, Vector<double>, Matrix2D<double>, double>(dRho, ddRho, normalVector, metricTensor, detm);
        }

        private Vector<double> CalclulateDeltaKsi (Vector<double> xUpdated, Dictionary<int, Vector<double>> dRho, Dictionary<int, Vector<double>> ddRho, Matrix2D<double> A, Matrix2D<double> metricTensor, double detm)
        {
            double f1 = dRho[1].DotProduct(A * xUpdated);
            double f2 = dRho[2].DotProduct(A * xUpdated);
            double e = ddRho[12].DotProduct(A * xUpdated);
            Vector<double> f = new Vector<double>(new double[] { f1, f2 });
            double[,] matrixM = new double[,] { { metricTensor[1, 1], e - metricTensor[0, 1] }, { e - metricTensor[1, 0], metricTensor[0, 0] } };
            Matrix2D<double> interMatrix = new Matrix2D<double>(matrixM);

            double coef = 1 / (detm - Math.Pow(e, 2) + 2 * e * metricTensor[0, 1]);
            Vector<double> deltaKsi = coef * (interMatrix * f);
            return deltaKsi;
        }

        private Vector<double> CPP(Vector<double> xUpdated)
        {
            double ksi1 = 0.0;
            double ksi2 = 0.0;
            Vector<double> ksiVector = new Vector<double>(new double[] { ksi1, ksi2 });
            for (int i = 0; i < 100; i++)
            {
                Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>> shapeFunctions = CalculateShapeFunctions(ksi1, ksi2);
                Dictionary<int, double> N = shapeFunctions.Item1;
                Dictionary<int, double> dN = shapeFunctions.Item2;
                Dictionary<int, double> ddN = shapeFunctions.Item3;
                Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>> positionMatrices = PositionMatrices(N, dN, ddN);
                IMatrix2D<double> A = positionMatrices.Item1;
                Dictionary<int, IMatrix2D<double>> dA = positionMatrices.Item2;
                Dictionary<int, IMatrix2D<double>> ddA = positionMatrices.Item3;
                Tuple<Dictionary<int, Vector<double>>, Dictionary<int, Vector<double>>, Vector<double>, Matrix2D<double>, double> surfaceProperties = SurfaceGeometry(xUpdated, dA, ddA);
                Dictionary<int, Vector<double>> dRho = surfaceProperties.Item1;
                Dictionary<int, Vector<double>> ddRho = surfaceProperties.Item2;
                Vector<double> n = surfaceProperties.Item3;
                Matrix2D<double> m = surfaceProperties.Item4;
                double detm = surfaceProperties.Item5;
                Vector<double> deltaKsi = CalclulateDeltaKsi(xUpdated, dRho, ddRho, (Matrix2D<double>)A, m, detm);
                if (deltaKsi.Norm<0.00001)
                {
                    break;
                }
                double[] ksi = ksiVector + deltaKsi;
                ksiVector = new Vector<double>(ksi); 
            }
            return ksiVector;
        }

        private double Penetration(Matrix2D<double> A, Vector<double> xUpdated, Vector<double> n)
        {
            double ksi3 = xUpdated.DotProduct(A * n);
            return ksi3;
        }

        private bool CheckContactStatus(Vector<double> xUpdated)
        {
            bool contactStatus;
            Vector<double> ksi = CPP(xUpdated);
            if (Math.Abs(ksi[0])>=1.05 && Math.Abs(ksi[1])>=1.05 )
            {
                Console.WriteLine("Projection point out of surface. No contact occurs");
                contactStatus = false;
                return contactStatus;
            }
            Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>> shapeFunctions = CalculateShapeFunctions(ksi[0], ksi[1]);
            Dictionary<int, double> N = shapeFunctions.Item1;
            Dictionary<int, double> dN = shapeFunctions.Item2;
            Dictionary<int, double> ddN = shapeFunctions.Item3;
            Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>> positionMatrices = PositionMatrices(N, dN, ddN);
            IMatrix2D<double> A = positionMatrices.Item1;
            Dictionary<int, IMatrix2D<double>> dA = positionMatrices.Item2;
            Dictionary<int, IMatrix2D<double>> ddA = positionMatrices.Item3;
            Tuple<Dictionary<int, Vector<double>>, Dictionary<int, Vector<double>>, Vector<double>, Matrix2D<double>, double> surfaceProperties = SurfaceGeometry(xUpdated, dA, ddA);
            Dictionary<int, Vector<double>> dRho = surfaceProperties.Item1;
            Dictionary<int, Vector<double>> ddRho = surfaceProperties.Item2;
            Vector<double> n = surfaceProperties.Item3;
            Matrix2D<double> m = surfaceProperties.Item4;
            double detm = surfaceProperties.Item5;
            double ksi3 = Penetration((Matrix2D<double>)A, xUpdated, n);
            if (ksi3 >= 0.0)
            {
                Console.WriteLine("No penetration. No contact occurs");
                contactStatus = false;
                return contactStatus;
            }
            contactStatus = true;
            return contactStatus;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element, Vector<double> xUpdated)
        {
            bool activeContact = CheckContactStatus(xUpdated);
            if (activeContact == false)
            {
                IMatrix2D<double> stiffnessMatrix = new Matrix2D<double>(new double[15, 15]);
            }
            Matrix2D<double> posA = (Matrix2D<double>)A;
            Matrix2D<double> mainPart = (posA.Transpose()*(normalVector ^ normalVector)*posA);
            return mainPart;
        }
    }
}
