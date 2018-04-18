﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Materials
{
    public class ElasticMaterial3D : IFiniteElementMaterial3D
    {
        //private readonly double[] strains = new double[6];
        private readonly double[] stresses = new double[6];
        private double[,] constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }
        private readonly double[] incrementalStrains = new double[6];
        private double[] stressesNew = new double[6];

        private double[,] GetConstitutiveMatrix()
        {
            double fE1 = YoungModulus / (double)(1 + PoissonRatio);
            double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
            double fE3 = fE1 + fE2;
            double fE4 = fE1 * 0.5;
            double[,] afE = new double[6, 6];
            afE[0, 0] = fE3;
            afE[0, 1] = fE2;
            afE[0, 2] = fE2;
            afE[1, 0] = fE2;
            afE[1, 1] = fE3;
            afE[1, 2] = fE2;
            afE[2, 0] = fE2;
            afE[2, 1] = fE2;
            afE[2, 2] = fE3;
            afE[3, 3] = fE4;
            afE[4, 4] = fE4;
            afE[5, 5] = fE4;

            //Vector<double> s = (new Matrix2D<double>(afE)) * (new Vector<double>(incrementalStrains));
            //s.Data.CopyTo(stresses, 0);

            return afE;
        }

        private void CalculateNextStressStrainPoint()
        {
            var stressesElastic = new double[6];
            for (int i = 0; i < 6; i++)
            {
                stressesElastic[i] = this.stresses[i];
                for (int j = 0; j < 6; j++)
                    stressesElastic[i] += this.constitutiveMatrix[i, j] * this.incrementalStrains[j];
            }

            this.stressesNew = stressesElastic;
        }

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { return 1; }
        }

        public bool Modified
        {
            get { return false; }
        }

        public void ResetModified()
        {
        }

        #endregion

        #region IFiniteElementMaterial3D Members

        public double[] Stresses { get { return stressesNew; } }

        public IMatrix2D<double> ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
                return new Matrix2D<double>(constitutiveMatrix);
            }
        }

        public void UpdateMaterial(double[] strainsIncrement)
        {
            //throw new NotImplementedException();
            Array.Copy(strainsIncrement, this.incrementalStrains, 6);
            constitutiveMatrix = GetConstitutiveMatrix();
            this.CalculateNextStressStrainPoint();

        }

        public void ClearState()
        {
            //Array.Clear(constitutiveMatrix, 0, constitutiveMatrix.Length);
            Array.Clear(incrementalStrains, 0, incrementalStrains.Length);
            Array.Clear(stresses, 0, stresses.Length);
            Array.Clear(stressesNew, 0, stressesNew.Length);
        }

        public void SaveState()
        {
            Array.Copy(this.stressesNew, this.stresses, 6);
        }

        public void ClearStresses()
        {
            Array.Clear(this.stresses, 0, 6);
            Array.Clear(this.stressesNew, 0, 6);
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new ElasticMaterial3D() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        #endregion

    }

}
