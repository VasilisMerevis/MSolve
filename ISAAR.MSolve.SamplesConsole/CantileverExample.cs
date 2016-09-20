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
    public static class CantileverExample
    {
        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 0.3, Y = 0 };
            Node node3 = new Node { ID = 3, X = 0.6, Y = 0 };
            Node node4 = new Node { ID = 4, X = 0.9, Y = 0 };
            Node node5 = new Node { ID = 5, X = 1.2, Y = 0 };
            Node node6 = new Node { ID = 6, X = 1.5, Y = 0 };
            Node node7 = new Node { ID = 7, X = 1.8, Y = 0 };
            Node node8 = new Node { ID = 8, X = 2.1, Y = 0 };
            Node node9 = new Node { ID = 9, X = 2.4, Y = 0 };
            Node node10 = new Node { ID = 10, X = 2.7, Y = 0 };
            Node node11 = new Node { ID = 11, X = 3, Y = 0 };

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

            return nodes;
        }
    }
}
