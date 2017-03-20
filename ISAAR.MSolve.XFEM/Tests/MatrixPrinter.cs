﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Tests
{
    static class MatrixPrinter
    {
        public static void PrintElementMatrix(int elementId, Matrix2D<double> k)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("K = \n" + k);
            Console.WriteLine("\n");
        }

        public static void PrintElementMatrices(int elementId, SymmetricMatrix2D<double> kss,
            Matrix2D<double> kes, SymmetricMatrix2D<double> kee)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("Kss = \n" + kss);
            Console.WriteLine("Kes = \n" + kes);
            Console.WriteLine("Kee = \n" + kee);
            Console.WriteLine("\n");
        }

        public static void PrintGlobalMatrix(Matrix2D<double> matrix)
        {
            Console.WriteLine("Global matrix:");
            Console.WriteLine("K = \n" + matrix);
            Console.WriteLine("\n");
        }
    }
}
