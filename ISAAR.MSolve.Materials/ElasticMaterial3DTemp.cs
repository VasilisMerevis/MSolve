﻿using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;

namespace ISAAR.MSolve.FEM.Materials
{
    public class ElasticMaterial3DTemp : IIsotropicContinuumMaterial3D
    {
        private readonly double[] strains = new double[6];
        private double[] incrementalStrains = new double[6];
        private double[] stresses = new double[6];
        private double[] stressesNew = new double[6];
        private double[,] constitutiveMatrix = null;
        public double YoungModulus { get; set; }
        public double PoissonRatio { get; set; }
        public double[] Coordinates { get; set; }

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

            //DIORTHOSI A.6
            //Vector s = (new Matrix2D(afE)) * (new Vector(strains));
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

            //this.stressesNew = stressesElastic; //DIORTHOSI A.4
            //Array.Copy(stressesElastic, this.stressesNew, 6);
            stressesNew = new double[6] { stressesElastic[0], stressesElastic[1], stressesElastic[2],
                stressesElastic[3], stressesElastic[4], stressesElastic[5] };

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

        public StressStrainVectorContinuum3D Stresses { get { return new StressStrainVectorContinuum3D(stressesNew); } } // DIORTHOSI A.3

        public ElasticityTensorContinuum3D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new StressStrainVectorContinuum3D(new double[6]));
                return new ElasticityTensorContinuum3D(constitutiveMatrix);
            }
        }

        public void UpdateMaterial(StressStrainVectorContinuum3D strainsIncrement)
        {
            //Array.Copy(strainsIncrement.Data, this.incrementalStrains, 6); //DIORHOSI A.5
            this.incrementalStrains = new double[6] { strainsIncrement.Data[0], strainsIncrement.Data[1], strainsIncrement.Data[2],
                strainsIncrement.Data[3], strainsIncrement.Data[4], strainsIncrement.Data[5] };

            constitutiveMatrix = GetConstitutiveMatrix();
            this.CalculateNextStressStrainPoint();
        }

        public void ClearState()
        {
            if (constitutiveMatrix != null) Array.Clear(constitutiveMatrix, 0, constitutiveMatrix.Length);
            //Array.Clear(incrementalStrains, 0, incrementalStrains.Length);
            //Array.Clear(stresses, 0, stresses.Length);
            //Array.Clear(stressesNew, 0, stressesNew.Length);
        }

        public void SaveState()
        {
            //Array.Copy(this.stressesNew, this.stresses, 6);
            this.stresses = new double[6] { stressesNew[0], stressesNew[1], stressesNew[2],
                stressesNew[3], stressesNew[4], stressesNew[5] };
        }

        public void ClearStresses()
        {
            //Array.Clear(this.stresses, 0, 6);
            //Array.Clear(this.stressesNew, 0, 6);
        }

        #endregion

        #region ICloneable Members
        object ICloneable.Clone() => Clone();

        public ElasticMaterial3DTemp Clone()
        {
            return new ElasticMaterial3DTemp() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
        }

        #endregion

    }
}
