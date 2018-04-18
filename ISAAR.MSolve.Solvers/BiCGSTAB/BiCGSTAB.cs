using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.BiCGSTAB;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;
using System.Diagnostics;
using System.IO;

namespace ISAAR.MSolve.Solvers.Skyline
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
                if (((SkylineMatrix2D<double>)subdomain.Matrix).IsFactorized) continue;

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
                
                double[] x = ((Vector<double>)subdomain.Solution).Data;
                Vector<double> xVector = new Vector<double>(x);

                Vector<double> pVector = new Vector<double>(new double[x.Length]);
                Vector<double> vVector = new Vector<double>(new double[x.Length]);

                Matrix2D<double> A = ((Matrix2D<double>)subdomain.Matrix);
                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();

                Vector<double> bVector = (Vector<double>)subdomain.RHS;
                double[] r = bVector-(A * xVector);

                Vector<double> rVector = new Vector<double>(r);
                double rho1 = rVector * rVector;

                double rho0 = 1.0;
                double w = 1.0;
                double a = 1.0;
                double b = rho1 * accuracyDigits / (rho0 * w);

                double[] inter1 = (pVector - w * vVector);
                Vector<double> inter1Vector = new Vector<double>(inter1);

                Vector<double> inter2Vector = b * inter1Vector;
                double[] inter3 = rVector + inter2Vector;
                pVector = new Vector<double>(inter3);

                vVector = A * pVector;
                a = rho1 / (rVector * vVector);
                Vector<double> sVector = new Vector<double>(rVector - a * vVector);
                Vector<double> tVector = A * sVector;
                w = (tVector * sVector) / (tVector * tVector);
                rho0 = rho1;
                rho1 = -w * (rVector * tVector);
                xVector = new Vector<double>(xVector + a * pVector);
                rVector = new Vector<double>(sVector - w * tVector);

                A.Solve(subdomain.RHS, x);
                stopWatch.Stop();
                //DestroyAccuracy(subdomain);

                x = new double[A.Rows];
                //LessenAccuracy(1e-7);
                A.Solve(subdomain.RHS, x);
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
