﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices
{
    /// <summary>
    /// Square non-symmetric invertible matrix. LU factorization needs pivoting.
    /// </summary>
    class SquareInvertible
    {
        public const int order = 10;

        public static double[,] matrix = new double[,] {
            { 1.0000,    2.6667,    7.0000,    5.0000,    2.5000,    9.0000,    6.0000,    2.2500,    4.0000,    3.0000 },
            { 0.0000,    0.3333,    1.0000,   -2.5000,    5.5000,   -7.7500,    2.0000,    5.7500,    3.0000,    2.0000 },
            { 2.0000,    2.0000,    4.0000,    9.0000,    1.7500,    5.0000,    3.5000,    2.0000,    6.0000,    8.0000 },
            { 1.0000,    1.0000,    3.0000,    5.0000,    1.7500,    6.0000,    4.0000,    7.0000,    6.0000,    9.0000 },
            { 5.0000,    0.6667,    4.5000,    5.0000,    0.2500,    2.3333,    1.5000,    1.2500,    9.0000,    0.7500 },
            { 0.3333,    2.0000,    4.0000,    5.0000,    5.0000,    9.0000,    2.5000,    4.0000,    1.0000,    3.0000 },
            { 0.2500,    6.0000,    8.0000,    0.7500,    1.0000,    3.0000,    8.0000,    1.2500,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 3.5000,    9.0000,    9.0000,    2.0000,    8.0000,    4.0000,    2.5000,    6.0000,    7.0000,    9.0000 },
            { 3.0000,    7.0000,    9.0000,    4.0000,    2.6667,    9.0000,    8.0000,    5.0000,    7.0000,    3.0000 }};

        public static double[] lhs = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064, 4.3251, 0.1117, 4.0487, 2.6311, 2.6269 };
        public static double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static double[,] lower = new double[,] {
            { 1.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.700000000000000, 1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.200000000000000, 0.296878936778343,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.200000000000000, 0.101561996458584, -0.113249698822564,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.000000000000000, 0.039058700551134,  0.176818919005663, -0.484906542960240,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.050000000000000, 0.699220466618463,  0.844451573716728, -0.331295118269131, -0.863097189214223,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.600000000000000, 0.773437271117538,  0.406895919175001,  0.052792355932608, -0.684415331679188,  0.241574580261808,  1.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.400000000000000, 0.203123992917168,  0.231873657358311,  0.942630231346478, -0.024512562835827,  0.346259134817979, -0.586219649578724,  1.000000000000000,  0.000000000000000, 0.000000000000000},
            { 0.200000000000000, 0.101561996458584,  0.345123356180875,  0.393366691999442,  0.151258061047848, -0.165741667413086,  0.355479089484414, -0.373132913343019,  1.000000000000000, 0.000000000000000},
            { 0.066660000000000, 0.229167553739405,  0.540742649585241,  0.391641616972103,  0.581445473311093, -0.628021544221750,  0.073364914263990,  0.198333706840520, -0.659496414881728, 1.000000000000000} };

        public static double[,] upper = new double[,] {
            { 5.000000000000000, 0.666700000000000, 4.500000000000000,  5.000000000000000, 0.250000000000000,  2.333300000000000, 1.500000000000000,  1.250000000000000,  9.000000000000000,  0.750000000000000 },
            { 0.000000000000000, 8.533310000000000, 5.850000000000000, -1.500000000000000, 7.825000000000000,  2.366690000000000, 1.450000000000000,  5.125000000000000,  0.699999999999999,  8.475000000000000 },
            { 0.000000000000000, 0.000000000000000, 4.363258219846695,  4.445318405167514, 0.126922319709468,  7.830719589116065, 5.269525541671404,  0.478495449010993,  1.992184744255161,  0.333951010803546 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  6.655773965243499, 0.169651292192538,  5.179800873632784, 2.149507285667150,  6.283684233638220, -0.645479075235204,  1.027081931408493 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 5.254188422424477, -6.715339873935628, 2.053923221210496,  8.512216710419459,  2.307405929749960,  2.107967404794224 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, -9.464114793276034, 4.946127816512377,  6.628553101741168,  0.155913636719474, -2.085750939237665 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000, 3.931767755773501,  3.984307335236524,  1.823617652005219, -2.248396632918166 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000, -5.326013241595878,  3.475940545000401,  4.388758952026945 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000,  4.020824604788080,  9.242293465315168 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000,  3.279189768619025 } };

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = MatrixMKL.CreateFromArray(matrix);
            var x = DenseVector.CreateFromArray(lhs);
            DenseVector b = A.MultiplyRight(x);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.InternalData);
        }

    }
}
