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

       

        #endregion
    }
}
