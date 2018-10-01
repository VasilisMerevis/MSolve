﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

//J_0inv_hexa and detJ_0 can only be calculated during initialization (at the first configuration) and then cached
namespace ISAAR.MSolve.FEM.Interpolation
{
    public class JacobianHexa8Reverse
    {
        public static (double[][,] J_0inv_hexa, double[] detJ_0) GetJ_0invHexaAndDetJ_0(IReadOnlyList<Matrix2D> ll1_hexa, IList<INode> elementNodes, int nGaussPoints)
        {
            double[][,] J_0b_hexa; // exoume tosa [,] osa einai kai ta gpoints
            double[][,] J_0_hexa;
            double[][,] J_0inv_hexa;
            double[] detJ_0; //osa kai ta gpoints

            double[][] ox_i;
            ox_i = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
            }
            J_0b_hexa = new double[nGaussPoints][,];
            J_0_hexa = new double[nGaussPoints][,];
            J_0inv_hexa = new double[nGaussPoints][,];
            detJ_0 = new double[nGaussPoints];

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                J_0b_hexa[gpoint] = new double[8, 3];
                J_0_hexa[gpoint] = new double[3, 3];
                J_0inv_hexa[gpoint] = new double[3, 3];

                //
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0b_hexa[gpoint][m, n] = ox_i[m][n];
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0_hexa[gpoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_0_hexa[gpoint][m, n] += ll1_hexa[gpoint][m, p] * J_0b_hexa[gpoint][p, n];
                        }
                    }
                }

                //
                double det1 = J_0_hexa[gpoint][0, 0] *
                         ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2]));
                double det2 = J_0_hexa[gpoint][0, 1] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]));
                double det3 = J_0_hexa[gpoint][0, 2] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1]));
                double jacobianDeterminant = det1 - det2 + det3;
                if (jacobianDeterminant < 0)
                {
                    throw new InvalidOperationException("The Jacobian Determinant is negative.");
                }
                detJ_0[gpoint] = jacobianDeterminant;

                //
                J_0inv_hexa[gpoint][0, 0] = ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2])) *
                                    (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 1] = ((J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 2] = ((J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 0] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 1] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 2] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 0] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 1] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 1]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 0])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 2] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 1]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 1])) *
                                        (1 / detJ_0[gpoint]);
            }


            return (J_0inv_hexa, detJ_0);
            

        }
        public static (double[][,] J_0inv_hexa, double[] detJ_0) GetJ_0invHexaAndDetJ_0(double[][,] ll1_hexa, IList<INode> elementNodes, int nGaussPoints)
        {
            double[][,] J_0b_hexa; // exoume tosa [,] osa einai kai ta gpoints
            double[][,] J_0_hexa;
            double[][,] J_0inv_hexa;
            double[] detJ_0; //osa kai ta gpoints

            double[][] ox_i;
            ox_i = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { elementNodes[j].X, elementNodes[j].Y, elementNodes[j].Z, };
            }
            J_0b_hexa = new double[nGaussPoints][,];
            J_0_hexa = new double[nGaussPoints][,];
            J_0inv_hexa = new double[nGaussPoints][,];
            detJ_0 = new double[nGaussPoints];

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                J_0b_hexa[gpoint] = new double[8, 3];
                J_0_hexa[gpoint] = new double[3, 3];
                J_0inv_hexa[gpoint] = new double[3, 3];

                //
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0b_hexa[gpoint][m, n] = ox_i[m][n];
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0_hexa[gpoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_0_hexa[gpoint][m, n] += ll1_hexa[gpoint][m, p] * J_0b_hexa[gpoint][p, n];
                        }
                    }
                }

                //
                double det1 = J_0_hexa[gpoint][0, 0] *
                         ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2]));
                double det2 = J_0_hexa[gpoint][0, 1] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]));
                double det3 = J_0_hexa[gpoint][0, 2] *
                              ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1]));
                double jacobianDeterminant = det1 - det2 + det3;
                if (jacobianDeterminant < 0)
                {
                    throw new InvalidOperationException("The Jacobian Determinant is negative.");
                }
                detJ_0[gpoint] = jacobianDeterminant;

                //
                J_0inv_hexa[gpoint][0, 0] = ((J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][1, 2])) *
                                    (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 1] = ((J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][0, 2] = ((J_0_hexa[gpoint][0, 1] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 1] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 0] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 2]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 1] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][2, 2]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][1, 2] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 2]) - (J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 2])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 0] = ((J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][2, 1]) - (J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][1, 1])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 1] = ((J_0_hexa[gpoint][2, 0] * J_0_hexa[gpoint][0, 1]) - (J_0_hexa[gpoint][2, 1] * J_0_hexa[gpoint][0, 0])) *
                                        (1 / detJ_0[gpoint]);
                J_0inv_hexa[gpoint][2, 2] = ((J_0_hexa[gpoint][0, 0] * J_0_hexa[gpoint][1, 1]) - (J_0_hexa[gpoint][1, 0] * J_0_hexa[gpoint][0, 1])) *
                                        (1 / detJ_0[gpoint]);
            }


            return (J_0inv_hexa, detJ_0);


        }

        public static double[][,] Get_J_1(int nGaussPoints, double[][] tx_i,IReadOnlyList<Matrix2D> ll1_hexa)
        {
            double[,] J_1b = new double[8, 3];
            double[][,] J_1 = new double[nGaussPoints][,];

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                J_1[npoint] = new double[3, 3];
            }

            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    //ll2[m, n] = tu_i[m][n];
                    J_1b[m, n] = tx_i[m][n];
                }
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_1[npoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_1[npoint][m, n] += ll1_hexa[npoint][m, p] * J_1b[p, n];
                        }
                    }
                }
            }

            return J_1;
        }
        public static double[][,] Get_J_1(int nGaussPoints, double[][] tx_i,double [][,] ll1_hexa)
        {
            double[,] J_1b = new double[8, 3];
            double[][,] J_1 = new double[nGaussPoints][,];

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                J_1[npoint] = new double[3, 3];
            }

            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    //ll2[m, n] = tu_i[m][n];
                    J_1b[m, n] = tx_i[m][n];
                }
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_1[npoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_1[npoint][m, n] += ll1_hexa[npoint][m, p] * J_1b[p, n];
                        }
                    }
                }
            }

            return J_1;
        }
    }
}
