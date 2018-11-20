using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.SamplesConsole
{
    public class TwoQuadsExample
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 2.0, Y = 0 };
            Node node3 = new Node { ID = 3, X = 2.0, Y = 2.0 };
            Node node4 = new Node { ID = 4, X = 2.0, Y = 4.0 };
            Node node5 = new Node { ID = 5, X = 0.0, Y = 4.0 };
            Node node6 = new Node { ID = 6, X = 0.0, Y = 2.0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);

            return nodes;
        }

        public static void Run()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 200.0e9;
            double poisson = 0.3;
            double load = 1.0e9;


            IList<Node> nodes = TwoQuadsExample.CreateNodes();

            Model twoBlocks = new Model();

            twoBlocks.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            for (int i = 0; i < nodes.Count; i++)
            {
                twoBlocks.NodesDictionary.Add(i + 1, nodes[i]);
            }

            twoBlocks.NodesDictionary[1].Constraints.Add(DOFType.X);
            twoBlocks.NodesDictionary[1].Constraints.Add(DOFType.Y);
            twoBlocks.NodesDictionary[2].Constraints.Add(DOFType.X);
            twoBlocks.NodesDictionary[2].Constraints.Add(DOFType.Y);

            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngMod,
                PoissonRatio = poisson
            };


            var element1 = new Element() { ID = 1, ElementType = new Quad4(material) { Thickness = 1.0 } };
            var element2 = new Element() { ID = 2, ElementType = new Quad4(material) { Thickness = 1.0 } };

            element1.AddNode(twoBlocks.NodesDictionary[1]);
            element1.AddNode(twoBlocks.NodesDictionary[2]);
            element1.AddNode(twoBlocks.NodesDictionary[3]);
            element1.AddNode(twoBlocks.NodesDictionary[6]);

            element2.AddNode(twoBlocks.NodesDictionary[6]);
            element2.AddNode(twoBlocks.NodesDictionary[3]);
            element2.AddNode(twoBlocks.NodesDictionary[4]);
            element2.AddNode(twoBlocks.NodesDictionary[5]);

            twoBlocks.ElementsDictionary.Add(element1.ID, element1);
            twoBlocks.ElementsDictionary.Add(element2.ID, element2);

            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);

            twoBlocks.Loads.Add(new Load() { Amount = load, Node = twoBlocks.NodesDictionary[4], DOF = DOFType.X });
            twoBlocks.Loads.Add(new Load() { Amount = -load, Node = twoBlocks.NodesDictionary[4], DOF = DOFType.Y });
            twoBlocks.Loads.Add(new Load() { Amount = -load, Node = twoBlocks.NodesDictionary[5], DOF = DOFType.Y });

            twoBlocks.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[0] = new SkylineLinearSystem(0, twoBlocks.SubdomainsDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            ProblemStructural provider = new ProblemStructural(twoBlocks, linearSystems);
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(twoBlocks.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(twoBlocks.Subdomains[0]) };
            int totalDOFs = twoBlocks.TotalDOFs;
            var linearSystemsArray = new[] { linearSystems[0] };
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters,
                subdomainMappers, provider, 10, totalDOFs);

            //var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            //linearSystems[0] = new SkylineLinearSystem(0, twoBlocks.SubdomainsDictionary[0].Forces);
            //SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            //ProblemStructural provider = new ProblemStructural(twoBlocks, linearSystems);

            //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);


            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
                twoBlocks.NodalDOFsDictionary[4][DOFType.Y],
                twoBlocks.NodalDOFsDictionary[5][DOFType.Y]
                });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Displacements of Node 3, along axes X, Y:");
            Console.WriteLine(childAnalyzer.Logs[0][0]);
        }
    }
}

