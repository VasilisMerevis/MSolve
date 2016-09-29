﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    public class Rod2D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[2] { DOFType.X, DOFType.Y };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly IFiniteElementMaterial material;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public double SectionArea { get; set; }

        public Rod2D(IFiniteElementMaterial material)
        {
            this.material = material;
        }

        public Rod2D(IFiniteElementMaterial material, IFiniteElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        #region IElementType Members

        public int ID
        {
            get { return 1; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.TwoD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double c2 = c * c;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
            double s2 = s * s;
            double cs = c * s;
            double E = (material as ElasticMaterial).YoungModulus;
            double A = SectionArea;
            
            return dofEnumerator.GetTransformedMatrix(
                new Matrix2D<double> (new double[,]
                {
                    {A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
                    {A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
                    {-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
                    {-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
                }));
        }

        public IMatrix2D<double> MassMatrix(Element element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            
            double totalMass = Density * SectionArea * L;

            return new Matrix2D<double>(new double[,]
            {
                { totalMass/2,0,0,0},
                {0,totalMass/2,0,0 },
                {0,0,totalMass/2,0 },
                {0,0,0,totalMass/2 }
            });
        }

        public IMatrix2D<double> DampingMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector<double> accelerations = new Vector<double>(4);
            IMatrix2D<double> massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            double[] forces = new double[4];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region IStructuralFiniteElement Members

        public IFiniteElementMaterial Material
        {
            get { return material; }
        }

        #endregion

        #region IFiniteElement Members


        public bool MaterialModified
        {
            get { return material.Modified; }
        }

        public void ResetMaterialModified()
        {
            material.ResetModified();
        }

        #endregion

        #region IFiniteElement Members

        public void ClearMaterialState()
        {
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
