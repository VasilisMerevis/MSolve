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
    public class ContactNtNTwoBLocks
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 1.0, Y = 0 };
            Node node3 = new Node { ID = 3, X = 1.0, Y = 1.0 };
            Node node4 = new Node { ID = 4, X = 0, Y = 1.0 };

            Node node5 = new Node { ID = 1, X = 0, Y = 1.01 };
            Node node6 = new Node { ID = 2, X = 1.0, Y = 1.01 };
            Node node7 = new Node { ID = 3, X = 1.0, Y = 2.01 };
            Node node8 = new Node { ID = 4, X = 0, Y = 2.01 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            nodes.Add(node5);
            nodes.Add(node6);
            nodes.Add(node7);
            nodes.Add(node8);

            return nodes;
        }

        public static void Run()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 100000.0;
            double poisson = 0.3;
            double loadY = 50000.0;
            

            IList<Node> nodes = ContactNtNTwoBLocks.CreateNodes();

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
            twoBlocks.NodesDictionary[7].Constraints.Add(DOFType.X);
            twoBlocks.NodesDictionary[8].Constraints.Add(DOFType.X);

            StressState2D stress = new StressState2D();
            ElasticMaterial2D mat = new ElasticMaterial2D(stress);


            var element1 = new Element() { ID = 1, ElementType = new Quad4(mat) };
            var element2 = new Element() { ID = 2, ElementType = new Quad4(mat) };
            
            var element3 = new Element() { ID = 3, ElementType = new Contact2DNtN(youngMod) };
            var element4 = new Element() { ID = 4, ElementType = new Contact2DNtN(youngMod) };

            element1.AddNode(twoBlocks.NodesDictionary[1]);
            element1.AddNode(twoBlocks.NodesDictionary[2]);
            element1.AddNode(twoBlocks.NodesDictionary[3]);
            element1.AddNode(twoBlocks.NodesDictionary[4]);

            element2.AddNode(twoBlocks.NodesDictionary[5]);
            element2.AddNode(twoBlocks.NodesDictionary[6]);
            element2.AddNode(twoBlocks.NodesDictionary[7]);
            element2.AddNode(twoBlocks.NodesDictionary[8]);

            element3.AddNode(twoBlocks.NodesDictionary[4]);
            element3.AddNode(twoBlocks.NodesDictionary[5]);

            element4.AddNode(twoBlocks.NodesDictionary[3]);
            element4.AddNode(twoBlocks.NodesDictionary[6]);

            twoBlocks.ElementsDictionary.Add(element1.ID, element1);
            twoBlocks.ElementsDictionary.Add(element2.ID, element2);
            twoBlocks.ElementsDictionary.Add(element3.ID, element3);
            twoBlocks.ElementsDictionary.Add(element4.ID, element4);

            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);
            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element3.ID, element3);
            twoBlocks.SubdomainsDictionary[0].ElementsDictionary.Add(element3.ID, element4);

            twoBlocks.Loads.Add(new Load() { Amount = loadY, Node = twoBlocks.NodesDictionary[7], DOF = DOFType.Y });
            twoBlocks.Loads.Add(new Load() { Amount = loadY, Node = twoBlocks.NodesDictionary[8], DOF = DOFType.Y });

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
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
                twoBlocks.NodalDOFsDictionary[3][DOFType.Y],
                twoBlocks.NodalDOFsDictionary[4][DOFType.Y]
                });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Displacements of Node 3, along axes X, Y:");
            Console.WriteLine(childAnalyzer.Logs[0][0]);
        }
    }
}
