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
    class CantileverExampleOneElement
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 0.3, Y = 0 };
            
            nodes.Add(node1);
            nodes.Add(node2);
            
            return nodes;
        }

        public static void Cantilever2DExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 200e9;
            double poisson = 0.3;
            double load = -2000000;

            ElasticMaterial material = new ElasticMaterial() { YoungModulus = youngMod, PoissonRatio = poisson };

            IList<Node> nodes = CantileverExampleOneElement.CreateNodes();

            Model cantiModel = new Model();

            cantiModel.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            for (int i = 0; i < nodes.Count; i++)
            {
                cantiModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.Y);
            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.RotZ);

            var element1 = new Element() { ID = 1, ElementType = new Beam2DNL(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };

            element1.AddNode(cantiModel.NodesDictionary[1]);
            element1.AddNode(cantiModel.NodesDictionary[2]);

            cantiModel.ElementsDictionary.Add(element1.ID, element1);

            cantiModel.SubdomainsDictionary[1].ElementsDictionary.Add(element1.ID, element1);

            cantiModel.Loads.Add(new Load() { Amount = load, Node = cantiModel.NodesDictionary[2], DOF = DOFType.Y });

            cantiModel.ConnectDataStructures();
            SolverSkyline solution = new SolverSkyline(cantiModel);

            ProblemStructural provider = new ProblemStructural(cantiModel, solution.SubdomainsDictionary);
            NonLinearAnalyzerNewtonRaphsonNew childAnalyzer = NonLinearAnalyzerNewtonRaphsonNew.NonLinearAnalyzerWithFixedLoadIncrements(solution, solution.SubdomainsDictionary, provider, 10, cantiModel.TotalDOFs);
            childAnalyzer.StepForMatrixRebuild = 1;
            //Analyzers.NewtonRaphsonNonLinearAnalyzer childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solution, solution.SubdomainsDictionary, provider, 10, cantiModel.TotalDOFs);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, solution.SubdomainsDictionary);

            //childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //    cantiModel.NodalDOFsDictionary[2][DOFType.X],
            //    cantiModel.NodalDOFsDictionary[2][DOFType.Y],
            //    cantiModel.NodalDOFsDictionary[2][DOFType.RotZ]});

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Console.WriteLine("Writing results for node 11");
            //Console.WriteLine("Dof and Values for Displacement Y");
            //Console.WriteLine(childAnalyzer.Logs[1][0]);
        }

    }
}
