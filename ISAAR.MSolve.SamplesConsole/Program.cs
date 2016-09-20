using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;






using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Elements;


namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private static void SolveBuildingInNoSoilSmall()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, 1, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = DOFType.X });
            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);
            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        static void Main(string[] args)
        {
            //SolveBuildingInNoSoilSmall();

            VectorExtensions.AssignTotalAffinityCount();
            double youngMod = 200e9;
            double poisson = 0.3;
            double load = -2000000;
            double area = 0.01;

            ElasticMaterial material = new ElasticMaterial() { YoungModulus = youngMod, PoissonRatio = poisson };

            IList<Node> nodes = CantileverExample.CreateNodes();

            Model cantiModel = new Model();

            cantiModel.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            for (int i = 0; i < nodes.Count; i++)
            {
                cantiModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.Y);
            cantiModel.NodesDictionary[1].Constraints.Add(DOFType.RotZ);

            var element1 = new Element() { ID = 1, ElementType = new Beam2D(material) { Density=1, SectionArea=0.01, MomentOfInertia = 8.333e-6} };
            var element2 = new Element() { ID = 2, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element3 = new Element() { ID = 3, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element4 = new Element() { ID = 4, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element5 = new Element() { ID = 5, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element6 = new Element() { ID = 6, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element7 = new Element() { ID = 7, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element8 = new Element() { ID = 8, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element9 = new Element() { ID = 9, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };
            var element10 = new Element() { ID = 10, ElementType = new Beam2D(material) { Density = 1, SectionArea = 0.01, MomentOfInertia = 8.333e-6 } };

            

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

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { cantiModel.NodalDOFsDictionary[11][DOFType.Y] });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Console.WriteLine("Writing results for node 11");
            Console.WriteLine("Dof and Values for Displacement Y");
            Console.WriteLine(childAnalyzer.Logs[1][0]);
        }
        
    }
}
