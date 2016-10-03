using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    /// <summary>
    /// Class for the Newton Rapshon Analyzer.
    /// </summary>
    /// <author>Theofilos Manitaras</author>
    public class NonLinearAnalyzerNewtonRaphsonNew : IAnalyzer
    {
        private static readonly int MAX_ITERATIOMS = 100;

        private static readonly double MAX_RESIDUAL_NORM_ALLOWED = 1.0e20;

        public int StepForMatrixRebuild { get; set; } = 1;

        public double Tolerance { get; set; } = 1.0e-5;

        private INonLinearParentAnalyzer parentAnalyzer;

        private Dictionary<int, Vector<double>> displacementMap = new Dictionary<int, Vector<double>>();

        private Vector<double> globalRightHandSide;

        private Dictionary<int, Vector<double>> incrementalDisplacementMap = new Dictionary<int, Vector<double>>();

        private double incrementalLoadFactor;

        private int increments;

        private Dictionary<int, Vector<double>> iterativeDisplacementMap = new Dictionary<int, Vector<double>>();

        private double loadFactor;

        private IList<double> loadIncrementsList;

        private INonLinearProvider provider;

        private double rightHandSideNorm;

        private Dictionary<int, Vector<double>> rightHandSides = new Dictionary<int, Vector<double>>();

        private ISolver solver;

        private Dictionary<int, ISolverSubdomain> subdomains;

        public IAnalyzer ParentAnalyzer
        {
            get { return (IAnalyzer)parentAnalyzer; }
            set { parentAnalyzer = (INonLinearParentAnalyzer)value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public Dictionary<int, IAnalyzerLog[]> Logs
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Factory method to create an instance of <see cref="NonLinearAnalyzerNewtonRaphsonNew"/> with
        /// a fixed incremental load factor.
        /// </summary>
        /// <param name="solver">The linear equation system solver.</param>
        /// <param name="subdomains">The subdomain dictionary.</param>
        /// <param name="provider">The provider.</param>
        /// <param name="increments">The number of load increments.</param>
        /// <param name="freedomDegreeCount"></param>
        /// <returns>The newly created instance.</returns>
        public static NonLinearAnalyzerNewtonRaphsonNew NonLinearAnalyzerWithFixedLoadIncrements(ISolver solver,
            Dictionary<int, ISolverSubdomain> subdomains, INonLinearProvider provider, int increments,
            int freedomDegreeCount)
        {
            IList<double> loadIncrements = new List<double>(increments);
            double loadIncrementFactor = 1.0 / increments;
            for (int i = 0; i < increments; i++)
            {
                loadIncrements.Add(loadIncrementFactor);
            }

            return new NonLinearAnalyzerNewtonRaphsonNew(solver, subdomains, provider, increments, freedomDegreeCount, loadIncrements);
        }

        /// <summary>
        /// Factory method to create an instance of <see cref="NonLinearAnalyzerNewtonRaphsonNew"/> with
        /// a given list of incremental load factors.
        /// </summary>
        /// <param name="solver">The linear equation system solver.</param>
        /// <param name="subdomains">The subdomain dictionary.</param>
        /// <param name="provider">The provider.</param>
        /// <param name="freedomDegreeCount"></param>
        /// <param name="loadIncrements">The list of incremental load factors.</param>
        /// <returns>The newly created instance.</returns>
        public static NonLinearAnalyzerNewtonRaphsonNew nonLinearAnalyzerWithPrescribedIncrements(ISolver solver,
                Dictionary<int, ISolverSubdomain> subdomains, INonLinearProvider provider, int freedomDegreeCount,
                IList<double> loadIncrements)
        {
            int increments = loadIncrements.Count;
            return new NonLinearAnalyzerNewtonRaphsonNew(solver, subdomains, provider, increments, freedomDegreeCount, loadIncrements);
        }

        private NonLinearAnalyzerNewtonRaphsonNew(ISolver solver, Dictionary<int, ISolverSubdomain> subdomains,
            INonLinearProvider provider, int increments, int freedomDegreeCount,
            IList<double> loadIncrements)
        {
            this.solver = solver;
            this.subdomains = subdomains;
            this.provider = provider;
            this.increments = increments;
            this.globalRightHandSide = new Vector<double>(freedomDegreeCount);
            this.loadIncrementsList = loadIncrements;
            foreach (var subdomain in subdomains.Values)
            {
                this.displacementMap[subdomain.ID] = new Vector<double>(subdomain.RHS.Length);
            }
        }

        public void Solve()
        {
            this.loadFactor = 0.0;
            this.initializeInternalVectors();

            for (int increment = 0; increment < this.increments; increment++)
            {
                this.incrementalLoadFactor = this.loadIncrementsList[increment];
                this.loadFactor += this.incrementalLoadFactor;
                this.clearIncrementalSolutionVector();
                this.UpdateRightHandSide();

                for (int step = 0; step < MAX_ITERATIOMS; step++)
                {

                    this.solver.Solve();
                    double errorNorm = this.CalculateInternalRightHandSideNorm() / (this.incrementalLoadFactor * this.rightHandSideNorm);

                    if ((Double.IsNaN(errorNorm)) || (errorNorm >= MAX_RESIDUAL_NORM_ALLOWED))
                    {
                        throw new ArithmeticException("The redisual norm exceeded the maximum allowed value.");
                    }
                    if (Math.Abs(errorNorm) < this.Tolerance)
                    {
                        break;
                    }

                    this.SplitResidualForcesToSubdomains();
                    if (((step + 1) % this.StepForMatrixRebuild) == 0)
                    {
                        this.BuildMatrices();
                    }

                }
                this.SaveMaterialStateAndUpdateSolution();
                Console.WriteLine(this.displacementMap[1][28]);

            }
            this.copySolutionToSubdomains();
        }


        private double CalculateInternalRightHandSideNorm()
        {
            this.globalRightHandSide.Clear();
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                Vector<double> iterativeDisplacements = this.iterativeDisplacementMap[id];
                iterativeDisplacements.Clear();
                iterativeDisplacements.Add(subdomain.Solution as Vector<double>);

                Vector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                incrementalDisplacements.Add(iterativeDisplacements);

                Vector<double> totalDisplacements = new Vector<double>(incrementalDisplacements.Length);
                Vector<double> displacements = this.displacementMap[id];
                totalDisplacements.Add(displacements);
                totalDisplacements.Add(incrementalDisplacements);

                Vector<double> internalRightHandSide = subdomain.GetRHSFromSolution(displacements, incrementalDisplacements) as Vector<double>;
                this.provider.ProcessInternalRHS(subdomain, internalRightHandSide.Data, totalDisplacements.Data);

                if (this.parentAnalyzer != null)
                {
                    Vector<double> additionalRightHandside = new Vector<double>(this.parentAnalyzer.GetOtherRHSComponents(subdomain, totalDisplacements));
                    internalRightHandSide.Add(additionalRightHandside);
                }

                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                Vector<double> rightHandSide = this.rightHandSides[id];
                subdomainRightHandSide.Clear();
                subdomainRightHandSide.Multiply(this.loadFactor);
                subdomainRightHandSide.Subtract(internalRightHandSide);
                subdomain.SubdomainToGlobalVector(subdomainRightHandSide.Data, this.globalRightHandSide.Data);
            }
            return this.provider.RHSNorm(this.globalRightHandSide.Data);
        }


        private void clearIncrementalSolutionVector()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                this.incrementalDisplacementMap[subdomain.ID].Clear();
            }
        }

        private void ClearMaterialStresses()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                subdomain.ClearMaterialStresses();
            }
        }


        private void copySolutionToSubdomains()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                double[] displacementData = new double[this.displacementMap[id].Length];
                this.displacementMap[id].CopyTo(displacementData, 0);
                subdomain.Solution = new Vector<double>(displacementData);
            }
        }


        private void initializeInternalVectors()
        {
            this.globalRightHandSide.Clear();
            this.rightHandSides.Clear();
            this.incrementalDisplacementMap.Clear();
            this.iterativeDisplacementMap.Clear();

            foreach (var subdomain in this.subdomains.Values)
            {
                int vectorLength = subdomain.RHS.Length;
                int id = subdomain.ID;
                double[] vecRData = new double[vectorLength];
                subdomain.RHS.CopyTo(vecRData, 0);
                Vector<double> vecR = new Vector<double>(vecRData);
                this.rightHandSides[id] = vecR;
                this.incrementalDisplacementMap[id] = new Vector<double>(vectorLength);
                this.iterativeDisplacementMap[id] = new Vector<double>(vectorLength);
                subdomain.SubdomainToGlobalVector((subdomain.RHS as Vector<double>).Data, this.globalRightHandSide.Data);
            }

            this.rightHandSideNorm = this.provider.RHSNorm(this.globalRightHandSide.Data);
        }


        private void SaveMaterialStateAndUpdateSolution()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                subdomain.SaveMaterialState();
                int id = subdomain.ID;
                Vector<double> displacements = this.displacementMap[id];
                Vector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                displacements.Add(incrementalDisplacements);
            }
        }


        private void SplitResidualForcesToSubdomains()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                subdomainRightHandSide.Clear();
                subdomain.SplitGlobalVectorToSubdomain(this.globalRightHandSide.Data, subdomainRightHandSide.Data);
            }
        }


        private void UpdateRightHandSide()
        {
            foreach (var subdomain in this.subdomains.Values)
            {
                int id = subdomain.ID;
                IVector<double> displacements = this.displacementMap[id];
                IVector<double> incrementalDisplacements = this.incrementalDisplacementMap[id];
                Vector<double> internalRightHandSide = subdomain.GetRHSFromSolution(displacements, incrementalDisplacements) as Vector<double>;

                if (this.parentAnalyzer != null)
                {
                    Vector<double> additionalRightHandside = new Vector<double>(
                            this.parentAnalyzer.GetOtherRHSComponents(subdomain, displacements));
                    internalRightHandSide.Add(additionalRightHandside);
                }

                Vector<double> subdomainRightHandSide = subdomain.RHS as Vector<double>;
                Vector<double> rightHandSide = this.rightHandSides[id];
                subdomainRightHandSide.Clear();
                subdomainRightHandSide.Add(rightHandSide);
                subdomainRightHandSide.Multiply(this.loadFactor);
                subdomainRightHandSide.Subtract(internalRightHandSide);
                subdomain.SubdomainToGlobalVector(subdomainRightHandSide.Data, this.globalRightHandSide.Data);
            }
        }

        public void BuildMatrices()
        {
            if (this.parentAnalyzer == null)
                throw new InvalidOperationException("This Newton-Raphson non-linear analyzer has no parent.");

            this.parentAnalyzer.BuildMatrices();
            this.solver.Initialize();
        }

        public void Initialize()
        {
            this.solver.Initialize();
        }
    }
}