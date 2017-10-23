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

            ElasticMaterial material = new ElasticMaterial() { YoungModulus = youngMod, PoissonRatio = poisson };

            IList<Node> nodes = TwoBlocksInContactExample.CreateNodes();

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

            var element1 = new Element() { ID = 1, ElementType = new Hexa8(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element2 = new Element() { ID = 2, ElementType = new Hexa8(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element3 = new Element() { ID = 3, ElementType = new Hexa8(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            



            element1.AddNode(cantiModel.NodesDictionary[1]);
            element1.AddNode(cantiModel.NodesDictionary[2]);

            element2.AddNode(cantiModel.NodesDictionary[2]);
            element2.AddNode(cantiModel.NodesDictionary[3]);

            element3.AddNode(cantiModel.NodesDictionary[3]);
            element3.AddNode(cantiModel.NodesDictionary[4]);

            element4.AddNode(cantiModel.NodesDictionary[4]);
            element4.AddNode(cantiModel.NodesDictionary[5]);

            element5.AddNode(cantiModel.NodesDictionary[5]);
            element5.AddNode(cantiModel.NodesDictionary[6]);

            element6.AddNode(cantiModel.NodesDictionary[6]);
            element6.AddNode(cantiModel.NodesDictionary[7]);

            element7.AddNode(cantiModel.NodesDictionary[7]);
            element7.AddNode(cantiModel.NodesDictionary[8]);

            element8.AddNode(cantiModel.NodesDictionary[8]);
            element8.AddNode(cantiModel.NodesDictionary[9]);

            element9.AddNode(cantiModel.NodesDictionary[9]);
            element9.AddNode(cantiModel.NodesDictionary[10]);

            element10.AddNode(cantiModel.NodesDictionary[10]);
            element10.AddNode(cantiModel.NodesDictionary[11]);

            cantiModel.ElementsDictionary.Add(element1.ID, element1);
            cantiModel.ElementsDictionary.Add(element2.ID, element2);
            cantiModel.ElementsDictionary.Add(element3.ID, element3);
            cantiModel.ElementsDictionary.Add(element4.ID, element4);
            cantiModel.ElementsDictionary.Add(element5.ID, element5);
            cantiModel.ElementsDictionary.Add(element6.ID, element6);
            cantiModel.ElementsDictionary.Add(element7.ID, element7);
            cantiModel.ElementsDictionary.Add(element8.ID, element8);
            cantiModel.ElementsDictionary.Add(element9.ID, element9);
            cantiModel.ElementsDictionary.Add(element10.ID, element10);

            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element1.ID, element1);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element2.ID, element2);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element3.ID, element3);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element4.ID, element4);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element5.ID, element5);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element6.ID, element6);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element7.ID, element7);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element8.ID, element8);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element9.ID, element9);
            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element10.ID, element10);

            cantiModel.Loads.Add(new Load() { Amount = load, Node = cantiModel.NodesDictionary[11], DOF = DOFType.Y });

            cantiModel.ConnectDataStructures();
            SolverSkyline solution = new SolverSkyline(cantiModel);

            ProblemStructural provider = new ProblemStructural(cantiModel, solution.SubdomainsDictionary);

            Analyzers.LinearAnalyzer childAnalyzer = new LinearAnalyzer(solution, solution.SubdomainsDictionary);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, solution.SubdomainsDictionary);

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
                cantiModel.NodalDOFsDictionary[11][DOFType.X],
                cantiModel.NodalDOFsDictionary[11][DOFType.Y],
                cantiModel.NodalDOFsDictionary[11][DOFType.RotZ]});

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Writing results for node 11");
            Console.WriteLine("Dof and Values for Displacement Y");
            Console.WriteLine(childAnalyzer.Logs[1][0]);
        }
    }
}
