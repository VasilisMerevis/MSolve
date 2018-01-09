using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using System.IO;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using System.Diagnostics;

namespace ISAAR.MSolve.Analyzers
{
    class CentralDifferencesDynamicAnalyzer : IAnalyzer
    {
        private readonly double alpha, delta, timeStep, totalTime;
        
        private Dictionary<int, Vector<double>> rhs = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> uu = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> uum = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> uc = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> ucc = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> u = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> v = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> v1 = new Dictionary<int, Vector<double>>();
        private Dictionary<int, Vector<double>> v2 = new Dictionary<int, Vector<double>>();
        
        
        private readonly IDictionary<int, ISolverSubdomain> subdomains;
        private readonly IImplicitIntegrationProvider provider;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer;

        private void InitializeInternalVectors()
        {
            uu.Clear();
            uum.Clear();
            uc.Clear();
            ucc.Clear();
            u.Clear();
            v.Clear();
            v1.Clear();
            v2.Clear();
            rhs.Clear();

            foreach (ISolverSubdomain subdomain in subdomains.Values)
            {
                int dofs = subdomain.RHS.Length;
                uu.Add(subdomain.ID, new Vector<double>(dofs));
                uum.Add(subdomain.ID, new Vector<double>(dofs));
                uc.Add(subdomain.ID, new Vector<double>(dofs));
                ucc.Add(subdomain.ID, new Vector<double>(dofs));
                u.Add(subdomain.ID, new Vector<double>(dofs));
                v.Add(subdomain.ID, new Vector<double>(dofs));
                v1.Add(subdomain.ID, new Vector<double>(dofs));
                v2.Add(subdomain.ID, new Vector<double>(dofs));
                rhs.Add(subdomain.ID, new Vector<double>(dofs));

                // Account for initial conditions coming from a previous solution
                subdomain.Solution.CopyTo(v[subdomain.ID].Data, 0);
            }
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return null; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

       public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            InitializeInternalVectors();
            //InitializeMatrices();
            //InitializeRHSs();
            childAnalyzer.Initialize();
        }

        public void BuildMatrices()
        {

        }

        public void Solve()
        {

        }

        #endregion
    }
}
