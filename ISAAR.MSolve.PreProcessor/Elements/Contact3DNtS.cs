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
    class Contact3DNtS 
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

        public Contact3DNtS(IFiniteElementMaterial3D material)
        {
            
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

            double[,] dA12 = new double[,]
            {
                {-ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0, 0},
                {0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0, 0 },
                {0, 0, -ddN[112], 0, 0, -ddN[212], 0, 0, -ddN[312], 0, 0, -ddN[412], 0, 0, 0 }
            };


            Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>> positionMatrices = 
                new Tuple<IMatrix2D<double>, Dictionary<int, IMatrix2D<double>>, Dictionary<int, IMatrix2D<double>>>(AMatrix, dA, dA2Matrix);
            return positionMatrices;
        }

        private Tuple<Vector<double>, Vector<double>, Matrix2D<double>, Vector<double>> SurfaceGeometry(IMatrix2D<double> dA1Matrix, IMatrix2D<double> dA2Matrix, IVector<double> xUpdated)
        {
            
            Vector<double> dRho1 = -1.0 * ((Matrix2D<double>)dA1Matrix * (Vector<double>)xUpdated);
            Vector<double> dRho2 = -1.0 * ((Matrix2D<double>)dA2Matrix * (Vector<double>)xUpdated);

            Matrix2D<double> metricTensor = new Matrix2D<double>(new double[,]
                { { dRho1*dRho1, dRho1*dRho2}, {dRho2*dRho1, dRho2*dRho2 } }
                );

            double detm = metricTensor[0, 0] * metricTensor[1, 1] - metricTensor[1, 0] * metricTensor[0, 1];

            Matrix2D<double> inverseMetricTensor = new Matrix2D<double>(new double[,] {
                { metricTensor[1,1]/detm, -metricTensor[0,1]/detm},
                { -metricTensor[1,0]/detm, metricTensor[0,0]/detm }
            });

            Vector<double> normalVector = (1/Math.Sqrt(detm)) * Vector<double>.CrossProductInR3(dRho1, dRho2);

            return new Tuple<Vector<double>, Vector<double>, Matrix2D<double>, Vector<double>>(dRho1, dRho2, metricTensor, normalVector);
        }

        private Vector<double> CalclulateDeltaKsi (Vector<double> xUpdated, Vector<double> dRho1, Vector<double> dRho2, Matrix2D<double> A)
        {
            double f1 = dRho1.DotProduct(A * xUpdated);
            double f2 = dRho2.DotProduct(A * xUpdated);
            Vector<double> f = new Vector<double>(new double[] { f1, f2 });
        }

        private Vector<double> CPP(Vector<double> xUpdated)
        {
            double ksi1 = 0.0;
            double ksi2 = 0.0;
            Tuple<Dictionary<int, double>, Dictionary<int, double>, Dictionary<int, double>> shapeFunctions = CalculateShapeFunctions(ksi1, ksi2);
            Dictionary<int,double> N = shapeFunctions.Item1;
            Dictionary<int, double> dN = shapeFunctions.Item2;
            Dictionary<int, double> ddN = shapeFunctions.Item3;
            Tuple<IMatrix2D<double>, IMatrix2D<double>, IMatrix2D<double>> positionMatrices = PositionMatrices(N, dN);
            IMatrix2D<double> A = positionMatrices.Item1;
            IMatrix2D<double> dA1 = positionMatrices.Item2;
            IMatrix2D<double> dA2 = positionMatrices.Item3;
            Tuple<Vector<double>, Vector<double>, Matrix2D<double>, Vector<double>> surfaceProperties = SurfaceGeometry(dA1, dA2, xUpdated);
            Vector<double> dRho1 = surfaceProperties.Item1;
            Vector<double> dRho2 = surfaceProperties.Item2;
            Matrix2D<double> m = surfaceProperties.Item3;
            Vector<double> n = surfaceProperties.Item4;
            
            
        }
    }
}
