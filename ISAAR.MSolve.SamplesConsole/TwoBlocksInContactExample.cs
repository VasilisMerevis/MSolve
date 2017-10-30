using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class TwoBlocksInContactExample
    {
        public static IList<Node> CreateNodes()
        {
            //Block 1
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 1.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 1.0, Y = 0.0, Z = 1.0 };
            Node node4 = new Node { ID = 4, X = 0.0, Y = 0.0, Z = 1.0 };
            Node node5 = new Node { ID = 5, X = 0.0, Y = 1.0, Z = 0.0 };
            Node node6 = new Node { ID = 6, X = 1.0, Y = 1.0, Z = 0.0 };
            Node node7 = new Node { ID = 7, X = 1.0, Y = 1.0, Z = 1.0 };
            Node node8 = new Node { ID = 8, X = 0.0, Y = 1.0, Z = 1.0 };
            //Block 2
            Node node9 = new Node { ID = 9, X = 0.25, Y = 1.0, Z = 0.25 };
            Node node10 = new Node { ID = 10, X = 0.75, Y = 1.0, Z = 0.25 };
            Node node11 = new Node { ID = 11, X = 0.75, Y = 1.0, Z = 0.75 };
            Node node12 = new Node { ID = 12, X = 0.25, Y = 1.0, Z = 0.75 };
            Node node13 = new Node { ID = 13, X = 0.25, Y = 1.5, Z = 0.25 };
            Node node14 = new Node { ID = 14, X = 0.75, Y = 1.5, Z = 0.25 };
            Node node15 = new Node { ID = 15, X = 0.75, Y = 1.5, Z = 0.75 };
            Node node16 = new Node { ID = 16, X = 0.25, Y = 1.5, Z = 0.75 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);
            nodes.Add(node5);
            nodes.Add(node6);
            nodes.Add(node7);
            nodes.Add(node8);
            nodes.Add(node9);
            nodes.Add(node10);
            nodes.Add(node11);
            nodes.Add(node12);
            nodes.Add(node13);
            nodes.Add(node14);
            nodes.Add(node15);
            nodes.Add(node16);

            return nodes;
        }

        public static void Two3DBlocksInContact()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 200e9;
            double poisson = 0.3;
            double load = -2000000;

            ElasticMaterial3D material = new ElasticMaterial3D() { YoungModulus = youngMod, PoissonRatio = poisson };

            IList<Node> nodes = CreateNodes();

            Model blocksModel = new Model();

            blocksModel.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            for (int i = 0; i < nodes.Count; i++)
            {
                blocksModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            blocksModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            blocksModel.NodesDictionary[1].Constraints.Add(DOFType.Y);
            blocksModel.NodesDictionary[1].Constraints.Add(DOFType.Z);
            blocksModel.NodesDictionary[2].Constraints.Add(DOFType.X);
            blocksModel.NodesDictionary[2].Constraints.Add(DOFType.Y);
            blocksModel.NodesDictionary[2].Constraints.Add(DOFType.Z);
            blocksModel.NodesDictionary[3].Constraints.Add(DOFType.X);
            blocksModel.NodesDictionary[3].Constraints.Add(DOFType.Y);
            blocksModel.NodesDictionary[3].Constraints.Add(DOFType.Z);
            blocksModel.NodesDictionary[4].Constraints.Add(DOFType.X);
            blocksModel.NodesDictionary[4].Constraints.Add(DOFType.Y);
            blocksModel.NodesDictionary[4].Constraints.Add(DOFType.Z);

            var element1 = new Element() { ID = 1, ElementType = new Hexa8(material) };
            var element2 = new Element() { ID = 2, ElementType = new Hexa8(material) };
            var element3 = new Element() { ID = 3, ElementType = new Contact3DNtS(material) };




            element1.AddNode(blocksModel.NodesDictionary[4]);
            element1.AddNode(blocksModel.NodesDictionary[3]);
            element1.AddNode(blocksModel.NodesDictionary[2]);
            element1.AddNode(blocksModel.NodesDictionary[1]);
            element1.AddNode(blocksModel.NodesDictionary[8]);
            element1.AddNode(blocksModel.NodesDictionary[7]);
            element1.AddNode(blocksModel.NodesDictionary[6]);
            element1.AddNode(blocksModel.NodesDictionary[5]);

            element2.AddNode(blocksModel.NodesDictionary[12]);
            element2.AddNode(blocksModel.NodesDictionary[11]);
            element2.AddNode(blocksModel.NodesDictionary[10]);
            element2.AddNode(blocksModel.NodesDictionary[9]);
            element2.AddNode(blocksModel.NodesDictionary[16]);
            element2.AddNode(blocksModel.NodesDictionary[15]);
            element2.AddNode(blocksModel.NodesDictionary[14]);
            element2.AddNode(blocksModel.NodesDictionary[13]);

            element3.AddNode(blocksModel.NodesDictionary[5]);
            element3.AddNode(blocksModel.NodesDictionary[6]);
            element3.AddNode(blocksModel.NodesDictionary[7]);
            element3.AddNode(blocksModel.NodesDictionary[8]);
            element3.AddNode(blocksModel.NodesDictionary[9]);



            blocksModel.ElementsDictionary.Add(element1.ID, element1);
            blocksModel.ElementsDictionary.Add(element2.ID, element2);
            blocksModel.ElementsDictionary.Add(element3.ID, element3);

            blocksModel.SubdomainsDictionary[1].ElementsDictionary.Add(element1.ID, element1);
            blocksModel.SubdomainsDictionary[1].ElementsDictionary.Add(element2.ID, element2);
            blocksModel.SubdomainsDictionary[1].ElementsDictionary.Add(element3.ID, element3);

            blocksModel.Loads.Add(new Load() { Amount = load, Node = blocksModel.NodesDictionary[13], DOF = DOFType.Y });
            blocksModel.Loads.Add(new Load() { Amount = load, Node = blocksModel.NodesDictionary[14], DOF = DOFType.Y });
            blocksModel.Loads.Add(new Load() { Amount = load, Node = blocksModel.NodesDictionary[15], DOF = DOFType.Y });
            blocksModel.Loads.Add(new Load() { Amount = load, Node = blocksModel.NodesDictionary[16], DOF = DOFType.Y });

            blocksModel.ConnectDataStructures();
            SolverSkyline linearSolution = new SolverSkyline(blocksModel);

            ProblemStructural provider = new ProblemStructural(blocksModel, linearSolution.SubdomainsDictionary);
            NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(linearSolution, linearSolution.SubdomainsDictionary, provider, 10, blocksModel.TotalDOFs);
            //Analyzers.LinearAnalyzer childAnalyzer = new LinearAnalyzer(solution, solution.SubdomainsDictionary);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSolution.SubdomainsDictionary);
            childAnalyzer.SetMaxIterations = 1000;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
                blocksModel.NodalDOFsDictionary[13][DOFType.X],
                blocksModel.NodalDOFsDictionary[13][DOFType.Y] });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Writing results for node 13");
            Console.WriteLine("Dof and Values for Displacement Y");
            Console.WriteLine(childAnalyzer.Logs[1][0]);
        }
    }
}
