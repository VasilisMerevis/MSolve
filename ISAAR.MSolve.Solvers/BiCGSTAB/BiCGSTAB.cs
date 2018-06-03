using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;
using System.Diagnostics;
using System.IO;

namespace ISAAR.MSolve.Solvers.BiCGSTAB
{
    public class BiCGSTAB : ISolver
    {
        private string stringFormat;
        private int accuracyDigits;
        private readonly Model model;
        private readonly Dictionary<int, ISolverSubdomain> subdomainsDictionary;
       

        public BiCGSTAB(Model model)
        {
            this.model = model;
            subdomainsDictionary = new Dictionary<int, ISolverSubdomain>(model.SubdomainsDictionary.Count);

            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                subdomainsDictionary.Add(subdomain.ID, new BiCGSTABSubdomain(subdomain));
        }

        public Dictionary<int, ISolverSubdomain> SubdomainsDictionary
        {
            get { return subdomainsDictionary; }
        }

        #region ISolver Members

        public void Initialize()
        {
            if (model.SubdomainsDictionary.Count != 1) throw new InvalidOperationException("Skyline solver operates on one subdomain only.");
            foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            {
                //if (((SkylineMatrix2D<double>)subdomain.Matrix).IsFactorized) continue;

                List<Vector<double>> zems = new List<Vector<double>>();
                List<int> zemColumns = new List<int>();
                Matrix2D<double> m = (Matrix2D<double>)subdomain.Matrix;

                Stopwatch stopWatch = new Stopwatch();
            }
        }

        public void Solve()
        {
            foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            {
                
                Vector<double> xVector = new Vector<double>(((Vector<double>)subdomain.Solution).Data);
                Vector<double> pVector = new Vector<double>(xVector.Data.Length);
                Vector<double> vVector = new Vector<double>(xVector.Data.Length);

                Matrix2D<double> A = ((Matrix2D<double>)subdomain.Matrix);
                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();

                Vector<double> bVector = (Vector<double>)subdomain.RHS;
                Vector<double> rVector = new Vector<double>(bVector - (A * xVector));
                Vector<double> r0hatVector = rVector;
                
                double rho0 = 1.0;
                double w = 1.0;
                double a = 1.0;
                double rho1 = r0hatVector * rVector;
                double b;
                double[] x;

                Vector<double> sVector;
                Vector<double> tVector;
                int iters = 100;
                double converged;
                for (int i = 0; i < iters; i++)
                {                   
                    b = (rho1 / rho0) * (a / w);
                    pVector = new Vector<double>(rVector + b * (new Vector<double>((pVector - w * vVector))));
                    vVector = A * pVector;
                    a = rho1 / (r0hatVector * vVector);
                    sVector = new Vector<double>(rVector - a * vVector);
                    tVector = A * sVector;
                    w = (tVector * sVector) / (tVector * tVector);
                    rho0 = rho1;
                    rho1 = -w * (r0hatVector * tVector);
                    xVector = new Vector<double>(xVector + a * pVector);
                    rVector = new Vector<double>(sVector - w * tVector);
                    converged = rVector.Norm;
                    if (i==iters | converged <  0.0001)
                    {
                        break;
                    }
                }
                x = xVector.Data;
                
                stopWatch.Stop();
               
                var xVec = new Vector<double>(x);
                var y = xVec.Norm;

                //StreamWriter sw = File.CreateText(@"c:\fbsub.txt");
                //sw.WriteLine(stopWatch.ElapsedMilliseconds.ToString());
                //sw.Close();
            }
        }

        #endregion
    }
}
