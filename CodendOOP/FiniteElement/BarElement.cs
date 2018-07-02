using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using Accord.Math;

namespace CodendOOP
{
    class BarElement
    {
        //======================
        // fields
        //======================

        public int ID;
        public Node n1;
        public Node n2;
        public double L0 { get; }

        public double[] B { get; private set; }
        public double L { get; private set; }

        public int Label { get; set; }
        public BarMaterial Material { get; set; }
        public Towing Towing { get; set; }

        public int[] globalDOF = new int[6];
        public double[] localF = new double[6];
        public double[,] localK = new double[6, 6];

        //======================
        // constructors
        //======================

        public BarElement(int ID, Node n1, Node n2, double L0)
        {
            this.ID = ID;
            this.n1 = n1;
            this.n2 = n2;
            this.L0 = L0;
            UpdateElement();
        }

        //=====================
        // methods
        //=====================

        public void PrintInfo()
        {
            Console.Write("{0,-10:D}", ID);
            Console.Write("{0,-10:D}{1,-10:D}{2,-10:D}", n1.ID, n2.ID);
            Console.Write("{0,-10:F2}", L0);
            Console.Write("{0,-10:F2}", L);
        }

        /*updating*/

        public void UpdateElement()
        {
            B = new double[3] {
            n2.X - n1.X,
            n2.Y - n1.Y,
            n2.Z - n1.Z
            };

            L = Sqrt(Pow(B[0], 2) + Pow(B[1], 2) + Pow(B[2], 2));
        }

        public void UpdateNodeID(List<Node> NodeList)
        {
            if (n1.CopyOf != -1)
            {
                n1 = NodeList.Find(node => node.ID == n1.CopyOf);
            }

            if (n2.CopyOf != -1)
            {
                n2 = NodeList.Find(node => node.ID == n2.CopyOf);
            }

            globalDOF = GetGlobalDOF();
        }

        public int[] GetGlobalDOF()
        {
            return new int[6] {n1.GetDof(0),
                               n1.GetDof(1),
                               n1.GetDof(2),
                               n2.GetDof(0),
                               n2.GetDof(1),
                               n2.GetDof(2)};
        }

        /*derrivatives*/

        public double DiffB(int w, int t)
        {
            if (w == t % 3 && t >= 0 && t < 6)
            {
                return Pow(-1, 1 + t / 3);
            }
            return 0;
        }

        public double DiffL(int t)
        {
            if (t >= 0 && t < 6)
            {
                return Pow(-1, 1 + t / 3) * B[t % 3] / L;
            }
            return 0;
        }

        public double Diff2L(int i, int j)
        {
            if (i >= 0 && i < 6)
            {
                return Pow(-1, 1 + i / 3) * DiffUnitB(i % 3, j);
            }
            return 0;
        }

        public double DiffUnitB(int w, int t)
        {
            if (t >= 0 && t < 6)
            {
                return (DiffB(w, t) * L - B[w] * DiffL(t)) / Pow(L, 2);
            }
            return 0;
        }
        
        public double Diff2UnitB(int w, int i, int j)
        {
            // w - component
            // i - derrivative with respect to vector x1 y1 y1 x2 y2 z2 x3 y3 z3
            // j - derrivative with respect to vector x1 y1 y1 x2 y2 z2 x3 y3 z3
            return ( ( DiffB(w,j) * DiffL(i) - DiffB(w, i) * DiffL(j) - B[w] * Diff2L(i, j) ) * Pow(L,2) - 2 * L * DiffL(i) * ( DiffB(w, j) * L - B[w] * DiffL(j) )  ) / Pow(L, 4);
        }

        private double DiffAlphaDrag(int t)
        {
            double dAlpha = 0;
            double cosa = Towing.Vector.Dot(B) / (Towing.Speed * L);
            double sina = Sqrt(1 - Pow(cosa, 2));

            for (int i = 0; i < 3; i++)
            {
                dAlpha = dAlpha - (Towing.Vector[i] / (Towing.Speed * L * sina)) * (DiffB(i, t) - B[i] * DiffL(t) / L);
            }

            return dAlpha;
        }

        private double DiffE(int w, int t)
        {
            return 2 * (B[0] * DiffB(0, t) + B[1] * DiffB(1, t) + B[2] * DiffB(2, t)) * Towing.Vector[w] -
                    (Towing.Vector[0] * DiffB(0, t) + Towing.Vector[1] * DiffB(1, t) + Towing.Vector[2] * DiffB(2, t)) * B[w] -
                    B.Dot(Towing.Vector) * DiffB(w, t);
        }



        public double DiffAlphaBending(int t, BarElement NextBar)
        {
            // when t 0:2 only thisBar derrivatives work
            // when t 6:8 only nextBar derrivatives work
            // when t 3:5 both thisBar and nextBar derrivatives work
            
            double cosa = B.Dot(NextBar.B) / (L * NextBar.L);
            double sina = Sqrt(1 - Pow(cosa, 2));
            if (sina != 0)
            {
            return (-1 / sina) * ((DiffUnitB(0,t) *  NextBar.B[0] + DiffUnitB(1, t) * NextBar.B[1] + DiffUnitB(2, t) * NextBar.B[2]) / NextBar.L +
                                  (NextBar.DiffUnitB(0, t - 3) * B[0] + NextBar.DiffUnitB(1, t - 3) * B[1] + NextBar.DiffUnitB(2, t - 3) * B[2]) / L);
            }
            return 0;
        }

        public double Diff2AlphaBending(int i, int j, BarElement NextBar)
        {
            double cosa = B.Dot(NextBar.B) / (L * NextBar.L);
            double sina = Sqrt(1 - Pow(cosa, 2));

            double Part1 = (-1 / sina);
            double dPart1 = DiffAlphaBending(i, NextBar) * cosa / Pow(sina, 2);
            double Part2 =
                (DiffUnitB(0, j) * NextBar.B[0] + DiffUnitB(1, j) * NextBar.B[1] + DiffUnitB(2, j) * NextBar.B[2]) / NextBar.L +
                (NextBar.DiffUnitB(0, j - 3) * B[0] + NextBar.DiffUnitB(1, j - 3) * B[1] + NextBar.DiffUnitB(2, j - 3) * B[2]) / L;
            double dPart2 =
            (Diff2UnitB(0, i, j) * NextBar.B[0] + Diff2UnitB(1, i, j) * NextBar.B[1] + Diff2UnitB(2, i, j) * NextBar.B[2]) / NextBar.L +
            (NextBar.Diff2UnitB(0, i - 3, j - 3) * B[0] + NextBar.Diff2UnitB(1, i - 3, j - 3) * B[1] + NextBar.Diff2UnitB(2, i - 3, j - 3) * B[2]) / L +
            DiffUnitB(0, i) * NextBar.DiffUnitB(0, j - 3) + DiffUnitB(1, i) * NextBar.DiffUnitB(1, j - 3) + DiffUnitB(2, i) * NextBar.DiffUnitB(2, j - 3) +
            NextBar.DiffUnitB(0, i - 3) * DiffUnitB(0, j) + NextBar.DiffUnitB(1, i - 3) * DiffUnitB(1, j) + NextBar.DiffUnitB(2, i - 3) * DiffUnitB(2, j); // correct
            
            return Part1 * dPart2 + Part2 * dPart1; 
        }

        public double DiffNormC(int t, BarElement NextBar)
        {
            double normC = Sqrt(Pow(B[0] + NextBar.B[0], 2) + Pow(B[1] + NextBar.B[1], 2) + Pow(B[2] + NextBar.B[2], 2));
            return (- 1 + t / 3) * (B[t % 3] + NextBar.B[t % 3]) / normC;
        }

        public double DiffCurv(int t, BarElement NextBar)
        {
            double normC = Sqrt(Pow(B[0] + NextBar.B[0], 2) + Pow(B[1] + NextBar.B[1], 2) + Pow(B[2] + NextBar.B[2], 2)); // never changes between 2 bars potential waste of memory

            double p = (L + NextBar.L + normC) / 2; // never changes between 2 bars potential waste of memory

            double S = Sqrt(p * (p - L) * (p - NextBar.L) * (p - normC));   // never changes between 2 bars potential waste of memory

            double DiffNormC = (-1 + t / 3) * (B[t % 3] + NextBar.B[t % 3]) / normC;

            double DiffLLC = DiffL(t) * NextBar.L * normC + L * NextBar.DiffL(t - 3) * normC + L * NextBar.L * DiffNormC;

            double DiffS = (DiffL(t) * L * Pow(NextBar.L,2) + NextBar.DiffL(t - 3) * Pow(L,2) * NextBar.L - DiffL(t) * Pow(L,3) +
                  NextBar.DiffL(t - 3) * NextBar.L * Pow(normC,2) + DiffNormC * Pow(NextBar.L,2) * normC - NextBar.DiffL(t - 3) * Pow(NextBar.L,3) +
                  DiffL(t) * L * Pow(normC,2) + DiffNormC * Pow(L,2) * normC - DiffNormC * Pow(normC,3)) / (8 * S);

            return 4 * (DiffS * L * NextBar.L * normC - DiffLLC * S) / Pow(L * NextBar.L * normC, 2);
        }


        public double CoeffPressureDrag(double d, double un)
        {
            if (un == 0)
            {
                return 0;
            }

            double rho = 1025;      // water density
            double mu = 0.00108;    // dynamic viscosity

            double Re = rho * d * un / mu;

            double s = -0.077215665 + Log(8 / Re);

            if (Re <= 1)
            {
                return (1 - 0.87 * Pow(s, -2)) * (8 * PI) / (Re * s);
            }

            if (Re <= 30)
            {
                return 1.45 + 8.55 * Pow(Re, -0.9);
            }
            else
            {
                return 1.1 + 4 * Pow(Re, -0.5);
            }
        }

        public double CoeffFrictionDrag(double d, double un)
        {
            if (un == 0)
            {
                return 0;
            }

            double rho = 1025;      // water density
            double mu = 0.00108;    // dynamic viscosity

            double Re = rho * d * un / mu;

            return PI * mu * (0.55 * Pow(Re, 0.5) + 0.084 * Pow(Re, (2.0 / 3.0)));
        }



        public double[] GetWeightForces()
        {
            double rhoWater = 1025;
            double[] F = new double[6];
            F[2] = F[5] = - PI * 0.25 * Pow(Material.Diameter, 2) * L * (Material.Density - rhoWater) * 9.81 / 2;
            return F;
        }

        public double[] GetTensionForces()
        {
            double EA = Material.EA;

            if (L <= L0)
            {
                EA = 0.0125 * Material.EA;
            }

            double magF = EA * (L - L0) / L0;

            return new double[6]
                {
                    magF * B[0] / L, // if tension (magF > 0) force is in the direction B vector
                    magF * B[1] / L,
                    magF * B[2] / L,
                  - magF * B[0] / L,
                  - magF * B[1] / L,
                  - magF * B[2] / L
                };
        }

        public double[] GetDragForces()
        {
            double rhoWater = 1025;      // water density
            double[] F = new double[6];

            // attack angle on the bar
            double cosa = Towing.Vector.Dot(B) / (Towing.Speed * L);
            double sina = Sqrt(1 - Pow(cosa, 2));
            double Cd = CoeffPressureDrag(Material.Diameter, Towing.Speed * Abs(sina));
            double Cf = CoeffFrictionDrag(Material.Diameter, Towing.Speed * Abs(sina));

            // normal force magnitue
            double magFn = 0.5 * rhoWater * Cd * Material.Diameter * L0 * Pow(Towing.Speed * sina, 2);
            // tangential force magintudes
            double magFt = 0.5 * rhoWater * Cf * Cd * Material.Diameter * L0 * Pow(Towing.Speed * cosa, 2);

            // current vector perpendicular to U and V twines
            double[] E = B.Cross(Towing.Vector.Cross(B));
            double normE = Sqrt(Pow(E[0], 2) + Pow(E[1], 2) + Pow(E[2], 2));

            for (int i = 0; i < 6; i++)
            {
                // normal force on the bar
                if (normE > 0)
                {
                    F[i] = F[i] + (magFn * E[i % 3] / normE) / 2;
                }
                // tangential force on the bar
                if (cosa != 0)
                {
                    F[i] = F[i] + (magFt * cosa * B[i % 3] / (Abs(cosa) * L)) / 2;
                }
            }
            return F;
        }

        public double[] GetBendingForces(BarElement NextBar)
        {
            double BdotB = B.Dot(NextBar.B);
            double cosa = BdotB / (L * NextBar.L);
            double sina = Sqrt(1 - Pow(cosa, 2));

            double[] C = new double[3] {
                    B[0] + NextBar.B[0],
                    B[1] + NextBar.B[1],
                    B[2] + NextBar.B[2] };

            double normC = Sqrt(Pow(C[0], 2) + Pow(C[1], 2) + Pow(C[2], 2));

            double p = (L + NextBar.L + normC) / 2;
            double S = Sqrt(p * (p - L) * (p - NextBar.L) * (p - normC));
            double R = 0.25 * (L * NextBar.L * normC) / S;

            double[] F = new double[9];

            if (sina != 0 && S != 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    F[i] = Material.EI * (B[i] * BdotB / Pow(L, 2) - NextBar.B[i]) / (L * NextBar.L * R * sina);
                    F[6 + i] = Material.EI * (-NextBar.B[i] * BdotB / Pow(NextBar.L, 2) + B[i]) / (L * NextBar.L * R * sina);
                    F[3 + i] = -F[6 + i] - F[i];
                }
            }
            return F;
        }

        public double[] GetTotalForces()
        {
            double[] totalF = new double[6];
            double[] F1 = new double[6];
            double[] F2 = new double[6];
            double[] F3 = new double[6];

            if (Material.Density > 0)
            {
                F1 = GetWeightForces();
            }

            if (Material.EA > 0)
            {
                F2 = GetTensionForces();
            }

            if (Towing.Speed > 0 && Towing.IncludeNettingDrag)
            {
                F3 = GetDragForces();
            }

            for (int i = 0; i < 6; i++)
            {
                totalF[i] = F1[i] + F2[i] + F3[i];
            }

            return totalF;
        }

        public void UpdateElementTotalForce()
        {
            Array.Clear(localF, 0, 6);
            localF = GetTotalForces();
        }



        public double[,] GetTensionStiffness()
        {
            double[,] K = new double[6, 6];

            double EA = Material.EA;

            if (L <= L0)
            {
                EA = 0.0125 * Material.EA;
            }

            double magF = EA * (L - L0) / L0;

            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    K[i, j] = -Pow(-1, i / 3) * ((EA / L0) * DiffL(j) * B[i % 3] / L + magF * DiffUnitB(i % 3, j)); // - correspoinds to the K = - J
                }
            }

            return K;
        }
      
        public double[,] GetDragStiffness()
        {
            double rhoWater = 1025;      // water density
            double[,] K = new double[6, 6];

            // attack angle on the bar
            double cosa = Towing.Vector.Dot(B) / (Towing.Speed * L);
            double sina = Sqrt(1 - Pow(cosa, 2));
            double Cd = CoeffPressureDrag(Material.Diameter, Towing.Speed * Abs(sina));
            double Cf = CoeffFrictionDrag(Material.Diameter, Towing.Speed * Abs(sina));

            // normal force magnitue
            double magFn = 0.5 * rhoWater * Cd * Material.Diameter * L0 * Pow(Towing.Speed * sina, 2);
            // tangential force magintudes
            double magFt = 0.5 * rhoWater * Cf * Cd * Material.Diameter * L0 * Pow(Towing.Speed * cosa, 2);

            // current vector perpendicular to U and V twines
            double[] E = B.Cross(Towing.Vector.Cross(B));
            double normE = Sqrt(Pow(E[0], 2) + Pow(E[1], 2) + Pow(E[2], 2));

            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    // normal force on U twines - tangent stiffness matrix
                    if (normE > 0 && sina != 0)
                    {
                        K[i, j] = K[i, j] + (
                                            rhoWater * Cd * Material.Diameter * L0 * Pow(Towing.Speed, 2) * cosa * sina * DiffAlphaDrag(j) * E[i % 3] / normE +
                                            (magFn / Pow(normE, 2)) * (DiffE(i % 3, j) * normE - (E[i % 3] / normE) * (E[0] * DiffE(0, j) + E[1] * DiffE(1, j) + E[2] * DiffE(2, j)))
                                            ) / -2;
                    }
                    // tangential force on U twines - tangent stiffness matrix
                    if (cosa != 0 && sina != 0)
                    {

                        K[i, j] = K[i, j] + (
                                             -Cf * rhoWater * Cd * Material.Diameter * L0 * Pow(Towing.Speed, 2) * cosa * sina * DiffAlphaDrag(j) * (cosa * B[i % 3] / (Abs(cosa) * L)) +
                                            (magFt / (L * Abs(cosa))) * (cosa * DiffB(i % 3, j) - sina * DiffAlphaDrag(j) * B[i % 3]) -
                                            (magFt * cosa * B[i % 3] / Pow(Abs(cosa) * L, 2)) * (Abs(cosa) * (B[0] * DiffB(0, j) + B[1] * DiffB(1, j) + B[2] * DiffB(2, j)) / L - cosa * sina * L * DiffAlphaDrag(j) / Abs(cosa))
                                             ) / -2;
                    }
                }
            }
            return K;
        }

        public double[,] GetBendingStiffness(BarElement NextBar)
        {
            double[,] K = new double[9, 9];

            double BdotB = B.Dot(NextBar.B);
            double cosa = BdotB / (L * NextBar.L);
            double sina = Sqrt(1 - Pow(cosa, 2));

            double[] C = new double[3] {
                    B[0] + NextBar.B[0],
                    B[1] + NextBar.B[1],
                    B[2] + NextBar.B[2] };

            double normC = Sqrt(Pow(C[0], 2) + Pow(C[1], 2) + Pow(C[2], 2));

            double p = (L + NextBar.L + normC) / 2;
            double S = Sqrt(p * (p - L) * (p - NextBar.L) * (p - normC));
            double R = 0.25 * (L * NextBar.L * normC) / S;

            if (sina != 0 && S != 0)
            {
                for (int i = 0; i < 9; i++)
                {
                    for (int j = 0; j < 9; j++)
                    {
                        K[i, j] = Material.EI * DiffCurv(j, NextBar) * DiffAlphaBending(i, NextBar) + (Material.EI / R) * Diff2AlphaBending(i, j, NextBar);
                    }
                }
            }

            return K;
        }

        public double[,] GetTotalStiffness()
        {
            double[,] totalK = new double[6,6];
            double[,] K1 = new double[6,6];
            double[,] K2 = new double[6,6];

            if (Material.EA > 0)
            {
                K1 = GetTensionStiffness();
            }

            if (Towing.Speed > 0 && Towing.IncludeNettingDrag)
            {
                K2 = GetDragStiffness();
            }

            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    totalK[i, j] = K1[i, j] + K2[i, j];
                }               
            }
            return totalK;
        }

        public void UpdateElementTotalStiffness()
        {
            Array.Clear(localK, 0, 36);
            localK = GetTotalStiffness();
        }
    }
}
