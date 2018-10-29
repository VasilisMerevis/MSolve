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

namespace ISAAR.MSolve.SamplesConsole
{
    public class ContactNtNLinearTruss
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 5.0, Y = 0 };
            Node node3 = new Node { ID = 3, X = 5.1, Y = 0 };
            Node node4 = new Node { ID = 4, X = 10.1, Y = 0 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            return nodes;
        }

        public static void Run()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 1000.0;
            double poisson = 0.3;
            double loadX = 50.0;
            double sectionArea = 1.0;

            IList<Node> nodes = ContactNtNLinearTruss.CreateNodes();

            Model trussModel = new Model();

            trussModel.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            for (int i = 0; i < nodes.Count; i++)
            {
                trussModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            trussModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[1].Constraints.Add(DOFType.Y);         
            trussModel.NodesDictionary[2].Constraints.Add(DOFType.Y);
            trussModel.NodesDictionary[3].Constraints.Add(DOFType.Y);
            trussModel.NodesDictionary[4].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[4].Constraints.Add(DOFType.Y);

            var element1 = new Element() { ID = 1, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };
            var element2 = new Element() { ID = 2, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };
            var element3 = new Element() { ID = 3, ElementType = new Contact2DNtN(youngMod) };

            element1.AddNode(trussModel.NodesDictionary[1]);
            element1.AddNode(trussModel.NodesDictionary[2]);

            element2.AddNode(trussModel.NodesDictionary[3]);
            element2.AddNode(trussModel.NodesDictionary[4]);

            element3.AddNode(trussModel.NodesDictionary[2]);
            element3.AddNode(trussModel.NodesDictionary[3]);

            trussModel.ElementsDictionary.Add(element1.ID, element1);
            trussModel.ElementsDictionary.Add(element2.ID, element2);
            trussModel.ElementsDictionary.Add(element3.ID, element3);

            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);
            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element3.ID, element3);

            trussModel.Loads.Add(new Load() { Amount = loadX, Node = trussModel.NodesDictionary[2], DOF = DOFType.X });

            trussModel.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[0] = new SkylineLinearSystem(0, trussModel.SubdomainsDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            ProblemStructural provider = new ProblemStructural(trussModel, linearSystems);
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(trussModel.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(trussModel.Subdomains[0]) };
            int totalDOFs = trussModel.TotalDOFs;
            var linearSystemsArray = new[] { linearSystems[0] };
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, 
                subdomainMappers, provider, 10, totalDOFs);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] {
                trussModel.NodalDOFsDictionary[2][DOFType.X],
                trussModel.NodalDOFsDictionary[3][DOFType.X]
                });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Displacements of Node 3, along axes X, Y:");
            Console.WriteLine(childAnalyzer.Logs[0][0]);
        }
    }
}
