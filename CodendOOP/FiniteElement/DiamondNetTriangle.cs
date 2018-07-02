using System;
using System.Collections.Generic;
using static System.Math;
using static System.Console;
using Accord.Math;

namespace CodendOOP
{
    class DiamondNetTriangle : TriangleElement
    {
        //=========================
        // variables
        //=========================

        public double normU;
        public double normV;

        //=========================
        // constructor
        //=========================

        public DiamondNetTriangle(int ID, Node n1, Node n2, Node n3, double[] uCoord, double[] vCoord) : base(ID, n1, n2, n3, uCoord, vCoord)
        {

        }

        //=========================
        // methods
        //=========================

        

        public override void GetUVd()
        {
            d = (UCoord[1] - UCoord[0]) * (VCoord[2] - VCoord[0]) -
                (UCoord[2] - UCoord[0]) * (VCoord[1] - VCoord[0]);

            U[0] = (n2.X - n1.X) * (VCoord[2] - VCoord[0]) / d -
                   (n3.X - n1.X) * (VCoord[1] - VCoord[0]) / d;

            U[1] = (n2.Y - n1.Y) * (VCoord[2] - VCoord[0]) / d -
                   (n3.Y - n1.Y) * (VCoord[1] - VCoord[0]) / d;

            U[2] = (n2.Z - n1.Z) * (VCoord[2] - VCoord[0]) / d -
                   (n3.Z - n1.Z) * (VCoord[1] - VCoord[0]) / d;

            normU = Sqrt(Pow(U[0], 2) + Pow(U[1], 2) + Pow(U[2], 2));

            V[0] = (n3.X - n1.X) * (UCoord[1] - UCoord[0]) / d -
                   (n2.X - n1.X) * (UCoord[2] - UCoord[0]) / d;

            V[1] = (n3.Y - n1.Y) * (UCoord[1] - UCoord[0]) / d -
                   (n2.Y - n1.Y) * (UCoord[2] - UCoord[0]) / d;

            V[2] = (n3.Z - n1.Z) * (UCoord[1] - UCoord[0]) / d -
                   (n2.Z - n1.Z) * (UCoord[2] - UCoord[0]) / d;

            normV = Sqrt(Pow(V[0], 2) + Pow(V[1], 2) + Pow(V[2], 2));
        }

        /* derrivatives */

        private double DiffU(int w, int t)
        {
            /* derrivative of Uw with respect to t
             * where w = {x y z} and t = {x1 y1 z1 x2 y2 z2 x3 y3 z3} */

            // for most of the w and t combinations dU = 0
            double outValue = 0;

            // dUx/dx1 or dUy/dy1 or dUz/dz1
            if ((w == 0 && t == 0) || (w == 1 && t == 1) || (w == 2 && t == 2))
            {
                outValue = (VCoord[1] - VCoord[2]) / d;
            }

            // dUx/dx2 or dUy/dy2 or dUz/dz2
            if ((w == 0 && t == 3) || (w == 1 && t == 4) || (w == 2 && t == 5))
            {
                outValue = (VCoord[2] - VCoord[0]) / d;
            }

            // dUx/dx3 or dUy/dy3 or dUz/dz3
            if ((w == 0 && t == 6) || (w == 1 && t == 7) || (w == 2 && t == 8))
            {
                outValue = (VCoord[0] - VCoord[1]) / d;
            }

            return outValue;
        }

        private double DiffV(int w, int t)
        {
            /* derrivative of Vw with respect to t
             * where w = {x y z} and t = {x1 y1 z1 x2 y2 z2 x3 y3 z3} */

            // for most of the w and t combinations dV = 0
            double outValue = 0;

            // dVx/dx1 or dVy/dy1 or dVz/dz1
            if ((w == 0 && t == 0) || (w == 1 && t == 1) || (w == 2 && t == 2))
            {
                outValue = (UCoord[2] - UCoord[1]) / d;
            }

            // dVx/dx2 or dVy/dy2 or dVz/dz2
            if ((w == 0 && t == 3) || (w == 1 && t == 4) || (w == 2 && t == 5))
            {
                outValue = (UCoord[0] - UCoord[2]) / d;
            }

            // dVx/dx3 or dVy/dy3 or dVz/dz3
            if ((w == 0 && t == 6) || (w == 1 && t == 7) || (w == 2 && t == 8))
            {
                outValue = (UCoord[1] - UCoord[0]) / d;
            }

            return outValue;
        }

        private double DiffNormU(int t)
        {
            /* derrivatives of normU with respect to t
             * where t = {x1 y1 z1 x2 y2 z2 x3 y3 z3} */
            return DiffU(t % 3, t) * U[t % 3] / normU;
        }

        private double DiffNormV(int t)
    {
        /* derrivatives of normV with respect to t
         * where t = {x1 y1 z1 x2 y2 z2 x3 y3 z3} */
        return DiffV(t % 3, t) * V[t % 3] / normV;
    }

        private double Diff2NormU(int i, int j)
    {
        double outValue = 0;
        if (i < 9 && j < 9)
        {
            outValue = DiffUnitU(0, i) * DiffU(0, j) + DiffUnitU(1, i) * DiffU(1, j) + DiffUnitU(2, i) * DiffU(2, j);
        }
        return outValue;
    }

        private double Diff2NormV(int i, int j)
        {
            double outValue = 0;
            if (i < 9 && j < 9)
            {
                outValue = DiffUnitV(0, i) * DiffV(0, j) + DiffUnitV(1, i) * DiffV(1, j) + DiffUnitV(2, i) * DiffV(2, j);
            }
            return outValue;
        }

        private double DiffUnitU(int w, int t)
        {
            double outValue = 0;
            if (w < 3 && t < 9)
            {
                outValue = (DiffU(w, t) * normU - U[w] * DiffNormU(t)) / Pow(normU, 2);
            }
            return outValue;
        }

        private double DiffUnitV(int w, int t)
        {
            double outValue = 0;
            if (w < 3 && t < 9)
            {
                outValue = (DiffV(w, t) * normV - V[w] * DiffNormV(t)) / Pow(normV, 2);
            }
            return outValue;
        }

        private double DiffAlphaDrag(int t)
        {
            double dAlpha = 0;
            double cosa = Towing.Vector.Dot(U) / (Towing.Speed * normU);
            double sina = Sqrt(1 - Pow(cosa, 2));

            if (sina != 0 && Towing.Speed > 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    dAlpha = dAlpha - (Towing.Vector[i] / (Towing.Speed * normU * sina)) * (DiffU(i, t) - U[i] * DiffNormU(t) / normU);
                }
            }

            return dAlpha;
        }

        private double DiffBetaDrag(int t)
        {
            double dBeta = 0;
            double cosb = Towing.Vector.Dot(V) / (Towing.Speed * normV);
            double sinb = Sqrt(1 - Pow(cosb, 2));

            if (sinb != 0 && Towing.Speed > 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    dBeta = dBeta - (Towing.Vector[i] / (Towing.Speed * normV * sinb)) * (DiffV(i, t) - V[i] * DiffNormV(t) / normV);
                }
            }

            return dBeta;
        }

        private double DiffEU(int w, int t)
        {
            double outValue = 0;

            if (w < 3 && t < 9)
            {
                outValue = 2 * (U[0] * DiffU(0, t) + U[1] * DiffU(1, t) + U[2] * DiffU(2, t)) * Towing.Vector[w] -
                            (Towing.Vector[0] * DiffU(0, t) + Towing.Vector[1] * DiffU(1, t) + Towing.Vector[2] * DiffU(2, t)) * U[w] -
                            U.Dot(Towing.Vector) * DiffU(w, t);
            }
            return outValue;
        }

        private double DiffEV(int w, int t)
        {
            double outValue = 0;

            if (w < 3 && t < 9)
            {
                outValue = 2 * (V[0] * DiffV(0, t) + V[1] * DiffV(1, t) + V[2] * DiffV(2, t)) * Towing.Vector[w] -
                            (Towing.Vector[0] * DiffV(0, t) + Towing.Vector[1] * DiffV(1, t) + Towing.Vector[2] * DiffV(2, t)) * V[w] -
                            V.Dot(Towing.Vector) * DiffV(w, t);
            }
            return outValue;
        }

        private double DiffAlphaOpen(int t)
        {
            double outValue = 0;
            double a = 0.5 * Acos(U.Dot(V) / (normU * normV));
            if (t < 9)
            {
                outValue = ((DiffUnitU(0, t) * V[0] + DiffUnitU(1, t) * V[1] + DiffUnitU(2, t) * V[2]) / normV +
                           (DiffUnitV(0, t) * U[0] + DiffUnitV(1, t) * U[1] + DiffUnitV(2, t) * U[2]) / normU) / (-2 * Sin(2 * a));
            }
            return outValue;
        }

        private double Diff2AlphaOpen(int i, int j)
        {
            double a = 0.5 * Acos(U.Dot(V) / (normU * normV));
            double outValue = 0;
            double A, dA, B, dB;

            if (i < 9 & j < 9)
            {
                A = DiffU(0, j) * V[0] + DiffU(1, j) * V[1] + DiffU(2, j) * V[2] - U.Dot(V) * DiffNormU(j) / normU +
                    DiffV(0, j) * U[0] + DiffV(1, j) * U[1] + DiffV(2, j) * U[2] - U.Dot(V) * DiffNormV(j) / normV;

                dA = DiffU(0, i) * DiffV(0, j) + DiffU(1, i) * DiffV(1, j) + DiffU(2, i) * DiffV(2, j) -
                        (DiffUnitU(0, i) * V[0] + DiffUnitU(1, i) * V[1] + DiffUnitU(2, i) * V[2]) * DiffNormU(j) -
                        (DiffV(0, i) * U[0] + DiffV(1, i) * U[1] + DiffV(2, i) * U[2]) * DiffNormU(j) / normU -
                        U.Dot(V) * Diff2NormU(i, j) / normU +
                        DiffV(0, i) * DiffU(0, j) + DiffV(1, i) * DiffU(1, j) + DiffV(2, i) * DiffU(2, j) -
                        (DiffUnitV(0, i) * U[0] + DiffUnitV(1, i) * U[1] + DiffUnitV(2, i) * U[2]) * DiffNormV(j) -
                        (DiffU(0, i) * V[0] + DiffU(1, i) * V[1] + DiffU(2, i) * V[2]) * DiffNormV(j) / normV -
                        U.Dot(V) * Diff2NormV(i, j) / normV;

                B = -2 * Sin(2 * a) * normU * normV;

                dB = -2 * Cos(2 * a) * 2 * DiffAlphaOpen(i) * normU * normV +
                     -2 * Sin(2 * a) * (DiffNormU(i) * normV + DiffNormV(i) * normU);

                outValue = (dA * B - dB * A) / Pow(B, 2);
            }
            return outValue;
        }

        /* force vectors */

        public override double[] GetNettingWeightForces()
        {
            double[] F = new double[9];
            F[2] = F[5] = F[8] = -d * PI * 0.25 * Pow(Material.TwineThickness, 2) * (Material.MeshSide / 2) * (Material.Density - 1025) * 9.81 / 3;
            return F;
        }

        public override double[] GetTwineTensionForces()
        {
            double[] F = new double[9];
            double l0 = Material.MeshSide / 2;
            // stiffnesses
            double EAU = Material.EA;
            double EAV = Material.EA;

            if (normU <= l0)
            {
                EAU = 0.0125 * Material.EA;
            }

            if (normV <= l0)
            {
                EAV = 0.0125 * Material.EA;
            }

            // F
            for (int i = 0; i < 3; i++)
            {
                F[0 * 3 + i] = (VCoord[2] - VCoord[1]) * (EAU * (normU - l0) / l0) * 0.5 * U[i] / normU +
                               (UCoord[1] - UCoord[2]) * (EAV * (normV - l0) / l0) * 0.5 * V[i] / normV;

                F[1 * 3 + i] = (VCoord[0] - VCoord[2]) * (EAU * (normU - l0) / l0) * 0.5 * U[i] / normU +
                               (UCoord[2] - UCoord[0]) * (EAV * (normV - l0) / l0) * 0.5 * V[i] / normV;

                F[2 * 3 + i] = (VCoord[1] - VCoord[0]) * (EAU * (normU - l0) / l0) * 0.5 * U[i] / normU +
                               (UCoord[0] - UCoord[1]) * (EAV * (normV - l0) / l0) * 0.5 * V[i] / normV;
            }

            return F;
        }

        public override double[] GetTwineDragForces()
        {
            double[] F = new double[9];

            double rhoWater = 1025;

            // attack angle on U twines
            double cosa = Towing.Vector.Dot(U) / (Towing.Speed * normU);
            double sina = Sqrt(1 - Pow(cosa, 2));
            double CdU = CoeffPressureDrag(Material.TwineThickness, Towing.Speed * Abs(sina));
            double CfU = CoeffFrictionDrag(Material.TwineThickness, Towing.Speed * Abs(sina));

            // attack angle on V twines
            double cosb = Towing.Vector.Dot(V) / (Towing.Speed * normV);
            double sinb = Sqrt(1 - Pow(cosb, 2));
            double CdV = CoeffPressureDrag(Material.TwineThickness, Towing.Speed * Abs(sinb));
            double CfV = CoeffFrictionDrag(Material.TwineThickness, Towing.Speed * Abs(sinb));

            // normal force magnitue
            double magFnU = 0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * sina, 2);
            double magFnV = 0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * sinb, 2);

            // tangential force magintudes
            double magFtU = CfU * 0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * cosa, 2);
            double magFtV = CfV * 0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * cosb, 2);

            // current vector perpendicular to U and V twines
            double[] EU = U.Cross(Towing.Vector.Cross(U));
            double[] EV = V.Cross(Towing.Vector.Cross(V));
            double normEU = Sqrt(Pow(EU[0], 2) + Pow(EU[1], 2) + Pow(EU[2], 2));
            double normEV = Sqrt(Pow(EV[0], 2) + Pow(EV[1], 2) + Pow(EV[2], 2));

            for (int i = 0; i < 9; i++)
            {
                // normal force on U twines
                if (normEU > 0)
                {
                    F[i] = F[i] + (magFnU * EU[i % 3] / normEU) / 3;
                }
                // normal force on V twines
                if (normEV > 0)
                {
                    F[i] = F[i] + (magFnV * EV[i % 3] / normEV) / 3;
                }
                // tangential force on U twines
                if (normU > 0 && cosa != 0)
                {
                    F[i] = F[i] + (magFtU * cosa * U[i % 3] / (Abs(cosa) * normU)) / 3;
                }
                // tangential force on V twines
                if (normV > 0 && cosb != 0)
                {
                    F[i] = F[i] + (magFtV * cosb * V[i % 3] / (Abs(cosb) * normV)) / 3;
                }
            }
            return F;
        }

        public override double[] GetCatchPressureForces()
        {
            double p = 0.5 * 1025 * 1.4 * Pow(Towing.Speed, 2);

            double[] Side12 = new double[] { n2.X - n1.X,
                                             n2.Y - n1.Y,
                                             n2.Z - n1.Z };

            double[] Side13 = new double[] { n3.X - n1.X,
                                             n3.Y - n1.Y,
                                             n3.Z - n1.Z };

            double[] crossProduct = Side12.Cross(Side13);

            return new double[]
                {
                            crossProduct[0] * p / 6,
                            crossProduct[1] * p / 6,
                            crossProduct[2] * p / 6,
                            crossProduct[0] * p / 6,
                            crossProduct[1] * p / 6,
                            crossProduct[2] * p / 6,
                            crossProduct[0] * p / 6,
                            crossProduct[1] * p / 6,
                            crossProduct[2] * p / 6
                };
        }

        public override double[] GetMeshOpeningForces()
        {
            double[] F = new double[9];

            double H = 0;
            if (orientation == 90)
            {
                H = Material.OpenningStifness;
            }
            else
            {
                H = - Material.OpenningStifness;
            }

            double a = 0.5 * Acos(U.Dot(V) / (normU * normV));

            for (int i = 0; i < 9; i++)
            {
                F[i] = - H * (a - Material.InitialOpeningAngle * PI / 180) * d * DiffAlphaOpen(i);
            }
            // contact between knots
            if (a <= Material.MinimumOpeningAngle * PI / 180)
            {
                for (int i = 0; i < 9; i++)
                {
                    F[i] = F[i] + Material.KnotContactStifness * (a - Material.MinimumOpeningAngle * PI / 180) * d * DiffAlphaOpen(i);
                }
            }
            return F;
        }

        public override double[] GetTwineBendingForces()
        {
            return new double[9];
        }

        public override double[] GetTotalForces()
        {
            double rhoWater = 1025;
            double[] totalF = new double[9];
            double[] F1 = new double[9];
            double[] F2 = new double[9];
            double[] F3 = new double[9];
            double[] F4 = new double[9];
            double[] F5 = new double[9];
            double[] F6 = new double[9];

            if (Material.Density != rhoWater)
            {
                F1 = GetNettingWeightForces();
            }

            if (Material.EA > 0)
            {
                F2 = GetTwineTensionForces();
            }

            if (Towing.Speed > 0) // check if current vector is nonzero
            {
                if (!HasCatch && Towing.IncludeNettingDrag)
                {
                    F3 = GetTwineDragForces();
                }
                
                if (HasCatch)
                {
                    F4 = GetCatchPressureForces();
                }
            }

            if (Material.OpenningStifness != 0 && orientation == 90)
            {
                F5 = GetMeshOpeningForces();
            }

            if (Material.EI > 0)
            {
                F6 = GetTwineBendingForces();
            }

            for (int i = 0; i < 9; i++)
            {
                totalF[i] = F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i];
            }
            //Write("\n");
            return totalF;
        }

        //public void UpdateElementTotalForce()
        //{
        //    Array.Clear(localF, 0, 9);
        //    localF = GetTotalForces();
        //}

        /* tangent stiffness matrices */

        public override double[,] GetTwineTensionStiffness()
        {
            double[,] K = new double[9, 9];
            double l0 = Material.MeshSide / 2;
            // stiffnesses
            double EAU = Material.EA;
            double EAV = Material.EA;

            if (normU <= l0)
            {
                EAU = 0.0125 * Material.EA;
            }

            if (normV <= l0)
            {
                EAV = 0.0125 * Material.EA;
            }

            // F
            for (int w = 0; w < 3; w++)
            {
                for (int t = 0; t < 9; t++)
                {
                    K[0 * 3 + w, t] = -0.5 * EAU * (VCoord[2] - VCoord[1]) * (DiffU(w, t) * (1 / l0 - 1 / normU) + DiffNormU(t) * U[w] / Pow(normU, 2))
                                     - 0.5 * EAV * (UCoord[1] - UCoord[2]) * (DiffV(w, t) * (1 / l0 - 1 / normV) + DiffNormV(t) * V[w] / Pow(normV, 2));

                    K[1 * 3 + w, t] = -0.5 * EAU * (VCoord[0] - VCoord[2]) * (DiffU(w, t) * (1 / l0 - 1 / normU) + DiffNormU(t) * U[w] / Pow(normU, 2))
                                     - 0.5 * EAV * (UCoord[2] - UCoord[0]) * (DiffV(w, t) * (1 / l0 - 1 / normV) + DiffNormV(t) * V[w] / Pow(normV, 2));

                    K[2 * 3 + w, t] = -0.5 * EAU * (VCoord[1] - VCoord[0]) * (DiffU(w, t) * (1 / l0 - 1 / normU) + DiffNormU(t) * U[w] / Pow(normU, 2))
                                     - 0.5 * EAV * (UCoord[0] - UCoord[1]) * (DiffV(w, t) * (1 / l0 - 1 / normV) + DiffNormV(t) * V[w] / Pow(normV, 2));
                }
            }

            return K;
        }

        public override double[,] GetTwineDragStiffness()
        {
            double[,] K = new double[9, 9];

            double rhoWater = 1025;

            // attack angle on U twines
            double cosa = Towing.Vector.Dot(U) / (Towing.Speed * normU);
            double sina = Sqrt(1 - Pow(cosa, 2));
            double CdU = CoeffPressureDrag(Material.TwineThickness, Towing.Speed * Abs(sina));
            double CfU = CoeffFrictionDrag(Material.TwineThickness, Towing.Speed * Abs(sina));

            // attack angle on V twines
            double cosb = Towing.Vector.Dot(V) / (Towing.Speed * normV);
            double sinb = Sqrt(1 - Pow(cosb, 2));
            double CdV = CoeffPressureDrag(Material.TwineThickness, Towing.Speed * Abs(sinb));
            double CfV = CoeffFrictionDrag(Material.TwineThickness, Towing.Speed * Abs(sinb));

            // normal force magnitue
            double magFnU = 0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * sina, 2);
            double magFnV = 0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * sinb, 2);

            // tangential force magintudes
            double magFtU = CfU * 0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * cosa, 2);
            double magFtV = CfV * 0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * 0.5 * d * Pow(Towing.Speed * cosb, 2);

            // current vector perpendicular to U and V twines
            double[] EU = U.Cross(Towing.Vector.Cross(U));
            double[] EV = V.Cross(Towing.Vector.Cross(V));
            double normEU = Sqrt(Pow(EU[0], 2) + Pow(EU[1], 2) + Pow(EU[2], 2));
            double normEV = Sqrt(Pow(EV[0], 2) + Pow(EV[1], 2) + Pow(EV[2], 2));

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    // normal force on U twines - tangent stiffness matrix
                    if (normEU > 0)
                    {
                        K[i, j] = K[i, j] + (
                                            0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * Pow(Towing.Speed, 2) * cosa * sina * d * DiffAlphaDrag(j) * EU[i % 3] / normEU +
                                            (magFnU / Pow(normEU, 2)) * (DiffEU(i % 3, j) * normEU - (EU[i % 3] / normEU) * (EU[0] * DiffEU(0, j) + EU[1] * DiffEU(1, j) + EU[2] * DiffEU(2, j)))
                                            ) / -3;
                    }
                    // normal force on V twines - tangent stiffness matrix
                    if (normEV > 0)
                    {
                        K[i, j] = K[i, j] + (
                                            0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * Pow(Towing.Speed, 2) * cosb * sinb * d * DiffBetaDrag(j) * EV[i % 3] / normEV +
                                            (magFnV / Pow(normEV, 2)) * (DiffEV(i % 3, j) * normEV - (EV[i % 3] / normEV) * (EV[0] * DiffEV(0, j) + EV[1] * DiffEV(1, j) + EV[2] * DiffEV(2, j)))
                                            ) / -3;
                    }
                    // tangential force on U twines - tangent stiffness matrix
                    if (normU > 0 && cosa != 0)
                    {

                        K[i, j] = K[i, j] + (
                                             -CfU * 0.5 * rhoWater * CdU * Material.TwineThickness * (Material.MeshSide / 2) * Pow(Towing.Speed, 2) * cosa * sina * d * DiffAlphaDrag(j) * (cosa * U[i % 3] / (Abs(cosa) * normU)) +
                                            (magFtU / (normU * Abs(cosa))) * (cosa * DiffU(i % 3, j) - sina * DiffAlphaDrag(j) * U[i % 3]) -
                                            (magFtU * cosa * U[i % 3] / Pow(Abs(cosa) * normU, 2)) * (Abs(cosa) * (U[0] * DiffU(0, j) + U[1] * DiffU(1, j) + U[2] * DiffU(2, j)) / normU - cosa * sina * normU * DiffAlphaDrag(j) / Abs(cosa))
                                             ) / -3;
                    }
                    // tangential force on V twines - tangent stiffness matrix
                    if (normV > 0 && cosb != 0)
                    {
                        K[i, j] = K[i, j] + (
                                            -CfV * 0.5 * rhoWater * CdV * Material.TwineThickness * (Material.MeshSide / 2) * Pow(Towing.Speed, 2) * cosb * sinb * d * DiffBetaDrag(j) * (cosb * V[i % 3] / (Abs(cosb) * normV)) +
                                            (magFtV / (normV * Abs(cosb))) * (cosb * DiffV(i % 3, j) - sinb * DiffBetaDrag(j) * V[i % 3]) -
                                            (magFtV * cosb * V[i % 3] / Pow(Abs(cosb) * normV, 2)) * (Abs(cosb) * (V[0] * DiffV(0, j) + V[1] * DiffV(1, j) + V[2] * DiffV(2, j)) / normV - cosb * sinb * normV * DiffBetaDrag(j) / Abs(cosb))
                                            ) / -3;
                    }
                }
            }

            //
            return K;
        }

        public override double[,] GetCatchPressureStiffness()
        {
            double p = 0.5 * 1025 * 1.4 * Pow(Towing.Speed, 2);

            double[] Side12 = new double[] { n2.X - n1.X,
                                             n2.Y - n1.Y,
                                             n2.Z - n1.Z };

            double[] Side13 = new double[] { n3.X - n1.X,
                                             n3.Y - n1.Y,
                                             n3.Z - n1.Z };

            double[] Side23 = new double[] { n3.X - n2.X,
                                             n3.Y - n2.Y,
                                             n3.Z - n2.Z };

            return new double[,]
            {
                { 0,          Side23[2] * p / 6, -Side23[1] * p / 6,    0        ,  -Side13[2] * p / 6, Side13[1] * p / 6,      0        , Side12[2] * p / 6,  -Side12[1] * p / 6 },
                {-Side23[2] * p / 6,  0        ,  Side23[0] * p / 6,   Side13[2] * p / 6,  0        ,  -Side13[0] * p / 6,     -Side12[2] * p / 6,  0        ,  Side12[0] * p / 6 },
                { Side23[1] * p / 6,  -Side23[0] * p / 6,  0        , -Side13[1] * p / 6,  Side13[0] * p / 6,  0        ,       Side12[1] * p / 6,  -Side12[0] * p / 6,  0        },
                { 0,          Side23[2] * p / 6, -Side23[1] * p / 6,    0        ,  -Side13[2] * p / 6, Side13[1] * p / 6,      0        , Side12[2] * p / 6,  -Side12[1] * p / 6 },
                {-Side23[2] * p / 6,  0        ,  Side23[0] * p / 6,   Side13[2] * p / 6,  0        ,  -Side13[0] * p / 6,     -Side12[2] * p / 6,  0        ,  Side12[0] * p / 6 },
                { Side23[1] * p / 6,  -Side23[0] * p / 6,  0        , -Side13[1] * p / 6,  Side13[0] * p / 6,  0        ,       Side12[1] * p / 6,  -Side12[0] * p / 6,  0        },
                { 0,          Side23[2] * p / 6, -Side23[1] * p / 6,    0        ,  -Side13[2] * p / 6, Side13[1] * p / 6,      0        , Side12[2] * p / 6,  -Side12[1] * p / 6 },
                {-Side23[2] * p / 6,  0        ,  Side23[0] * p / 6,   Side13[2] * p / 6,  0        ,  -Side13[0] * p / 6,     -Side12[2] * p / 6,  0        ,  Side12[0] * p / 6 },
                { Side23[1] * p / 6,  -Side23[0] * p / 6,  0        , -Side13[1] * p / 6,  Side13[0] * p / 6,  0        ,       Side12[1] * p / 6,  -Side12[0] * p / 6,  0        }
            };
        }

        public override double[,] GetMeshOpeningStiffness()
        {
            double H = 0;
            if (orientation == 90)
            {
                H = Material.OpenningStifness;
            }
            else
            {
                H =  -Material.OpenningStifness;
            }

            double[,] K = new double[9, 9];
            double a = 0.5 * Acos(U.Dot(V) / (normU * normV));

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    K[i, j] = H * d * DiffAlphaOpen(i) * DiffAlphaOpen(j) +
                              H * d * (a - Material.InitialOpeningAngle * PI /180) * Diff2AlphaOpen(i, j);
                }
            }
            return K;
        }

        public override double[,] GetTwineBendingStiffness()
        {
            return new double[9, 9];
        }

        public override double[,] GetTotalStiffness()
        {
            double[,] totalK = new double[9, 9];
            double[,] K2 = new double[9, 9];
            double[,] K3 = new double[9, 9];
            double[,] K4 = new double[9, 9];
            double[,] K5 = new double[9, 9];
            double[,] K6 = new double[9, 9];

            if (Material.EA > 0)
            {
                K2 = GetTwineTensionStiffness();
            }

            if (Towing.Speed != 0) // check if current vector is nonzero
            {
                if (!HasCatch && Towing.IncludeNettingDrag)
                {
                    K3 = GetTwineDragStiffness();
                }
                if (HasCatch)
                {
                    K4 = GetCatchPressureStiffness();
                }
            }

            if (Material.OpenningStifness != 0 && orientation == 90)
            {
                K5 = GetMeshOpeningStiffness();
            }

            if (Material.EI > 0)
            {
                K6 = GetTwineBendingStiffness();
            }

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    totalK[i, j] = K2[i, j] + K3[i, j] + K4[i, j] + K5[i, j] + K6[i, j];
                }
            }
            return totalK;
        }

        //public void UpdateElementTotalStiffness()
        //{
        //    Array.Clear(localK, 0, 81);
        //    localK = GetTotalStiffness();
        //}
    }
}
