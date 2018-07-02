using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.IO;
using CSparse;
using CSparse.Storage;
using CSparse.Double.Factorization;

namespace CodendOOP
{
    class Simulation
    {

        //=========================
        // variables
        //=========================

        public Codend Codend;
        public Catch Catch;
        public Towing Towing;

        public int totalDOFcount;
        public int freeDOFcount;
        public int fixedDOFcount;
        public int[] fixedDOF;
        public int[] freeDOF;

        public double[] globalF;
        public CoordinateStorage<double> globalK;
        public double residual;

        public SolverSettings SolverSettings { get; set; }
        public InputReader FilePaths { get; set; }

        //=========================
        // constructors
        //=========================

        public Simulation(Codend Codend, Catch Catch, Towing Towing)
        {
            this.Codend = Codend;
            this.Catch = Catch;
            this.Towing = Towing;

            SolverSettings = new SolverSettings(); // start with default settings

            Codend.SetCylInitialShape();
            Codend.RestrainEntrance();
            CalculateBoundaryDOF();
        }

        //=========================
        // methods
        //=========================

        private void InitializeFEM()
        {
            globalF = new double[totalDOFcount];
            int localKsize = 81;
            globalK = new CoordinateStorage<double>(totalDOFcount, totalDOFcount, localKsize * Codend.TriangleList.Count);
        }

        private void ApplyTowing()
        {
            for (int i = 0; i < Codend.TriangleList.Count; i++)
            {
                Codend.TriangleList[i].Towing = Towing;
            }

            if (Codend.BarList.Count != 0)
            {
                for (int i = 0; i < Codend.BarList.Count; i++)
                {
                    Codend.BarList[i].Towing = Towing;
                }
            }
        }

        private void CalculateBoundaryDOF()
        {
            int dofPerNode = 3;
            totalDOFcount = Codend.NodeList.Count * dofPerNode;
            
            if(totalDOFcount == 0)
            {
                throw new ArgumentException("Boundary conditions cannot be established for empty domain.");
            }

            int[] tempFixedDOF = new int[dofPerNode * Codend.RestraintList.Count];
            fixedDOFcount = 0;

            // make temporary fixedDOF array
            for (int i = 0; i < Codend.RestraintList.Count; i++)
            {
                // check if x is fixed
                if (Codend.RestraintList[i].xIsFixed)
                {
                    tempFixedDOF[3 * i] = 3 * Codend.RestraintList[i].nodeID;
                    fixedDOFcount++;
                }
                else
                {
                    tempFixedDOF[3 * i] = -1;
                }
                // check if y is fixed
                if (Codend.RestraintList[i].yIsFixed)
                {
                    tempFixedDOF[3 * i + 1] = 3 * Codend.RestraintList[i].nodeID + 1;
                    fixedDOFcount++;
                }
                else
                {
                    tempFixedDOF[3 * i + 1] = -1;
                }
                // check if z is fixed
                if (Codend.RestraintList[i].zIsFixed)
                {
                    tempFixedDOF[3 * i + 2] = 3 * Codend.RestraintList[i].nodeID + 2;
                    fixedDOFcount++;
                }
                else
                {
                    tempFixedDOF[3 * i + 2] = -1;
                }
            }

            // sort temporary fixedDOF array
            fixedDOF = new int[fixedDOFcount];
            int count = 0;
            for (int i = 0; i < tempFixedDOF.Length; i++)
            {
                if (tempFixedDOF[i] != -1)
                {
                    fixedDOF[count] = tempFixedDOF[i];
                    count++;
                }
            }

            // free DOF
            freeDOFcount = totalDOFcount - fixedDOFcount;
            freeDOF = new int[freeDOFcount];

            count = 0;

            for (int i = 0; i < totalDOFcount; i++)
            {
                if (IsFreeDOF(i))
                {
                    freeDOF[count] = i;
                    count++;
                }
            }
        }

        private bool IsFreeDOF(int dof)
        {
            bool isFree = true;
            for (int i = 0; i < fixedDOFcount; i++)
            {
                if (dof == fixedDOF[i])
                {
                    isFree = false;
                    break;
                }
            }
            return isFree;
        }

        private double UpdateForceResidual()
        {
            residual = 0;
            for (int i = 0; i < totalDOFcount; i++)
            {
                residual = residual + Math.Pow(globalF[i], 2);
            }
            return Math.Sqrt(residual);
        }

        private double DisplacmentNorm(double[] h)
        {
            double norm = 0;
            for (int i = 0; i < totalDOFcount; i++)
            {
                norm += Math.Pow(h[i], 2);                      // residue as sum of squares
            }
            return Math.Sqrt(norm);
        }

        private void AssembleForces()
        {
            InitializeFEM();

            /* TRIANGLES */
            Parallel.For(0, Codend.TriangleList.Count, elem =>
            {
                Codend.TriangleList[elem].UpdateElementTotalForce();
                Codend.TriangleList[elem].UpdateElementTotalStiffness();

                for (int row = 0; row < 9; row++)
                {
                    //if (fixedDOF.Any(x => x == Codend.TriangleList[elem].globalDOF[row]))
                    if (!IsFreeDOF(Codend.TriangleList[elem].globalDOF[row]))
                    {
                        continue;
                    }

                    lock (this)
                    {
                        globalF[Codend.TriangleList[elem].globalDOF[row]] += Codend.TriangleList[elem].localF[row];
                    }
                }
            });

            /* BARS */
            if (Codend.BarList.Count != 0)
            {
                Parallel.For(0, Codend.BarList.Count, elem =>
                {
                    Codend.BarList[elem].UpdateElementTotalForce();
                    Codend.BarList[elem].UpdateElementTotalStiffness();

                    for (int row = 0; row < 6; row++)
                    {
                        //if (fixedDOF.Any(x => x == Codend.BarList[elem].globalDOF[row]))
                        if (!IsFreeDOF(Codend.BarList[elem].globalDOF[row]))
                        {
                            continue;
                        }

                        lock (this)
                        {
                            globalF[Codend.BarList[elem].globalDOF[row]] += Codend.BarList[elem].localF[row];
                        }
                    }
                });
            }
        }

        private void AssembleForcesAndStiffness(double kDiag)
        {
            InitializeFEM();

            /* TRIANGLES */
            Parallel.For(0, Codend.TriangleList.Count, elem =>
            {
                Codend.TriangleList[elem].UpdateElementTotalForce();
                Codend.TriangleList[elem].UpdateElementTotalStiffness();

                for (int row = 0; row < 9; row++)
                {
                    //if (fixedDOF.Any(x => x == Codend.TriangleList[elem].globalDOF[row]))
                    if (!IsFreeDOF(Codend.TriangleList[elem].globalDOF[row]))
                    {
                        continue;
                    }

                    lock (this)
                    {
                        globalF[Codend.TriangleList[elem].globalDOF[row]] += Codend.TriangleList[elem].localF[row];
                    }

                    for (int col = 0; col < 9; col++)
                    {
                        //if (fixedDOF.Any(x => x == Codend.TriangleList[elem].globalDOF[col]))
                        if (!IsFreeDOF(Codend.TriangleList[elem].globalDOF[col]))
                        {
                            continue;
                        }

                        lock (this)
                        {
                            globalK.At(Codend.TriangleList[elem].globalDOF[row],
                                             Codend.TriangleList[elem].globalDOF[col],
                                             Codend.TriangleList[elem].localK[row, col]);
                        }

                    }
                }
            });

            /* BARS */
            if (Codend.BarList.Count != 0)
            {
                Parallel.For(0, Codend.BarList.Count, elem =>
                {
                    Codend.BarList[elem].UpdateElementTotalForce();
                    Codend.BarList[elem].UpdateElementTotalStiffness();

                    for (int row = 0; row < 6; row++)
                    {
                    //if (fixedDOF.Any(x => x == Codend.BarList[elem].globalDOF[row]))
                    if (!IsFreeDOF(Codend.BarList[elem].globalDOF[row]))
                        {
                            continue;
                        }

                        lock (this)
                        {
                            globalF[Codend.BarList[elem].globalDOF[row]] += Codend.BarList[elem].localF[row];
                        }

                        for (int col = 0; col < 6; col++)
                        {
                        //if (fixedDOF.Any(x => x == Codend.BarList[elem].globalDOF[col]))
                        if (!IsFreeDOF(Codend.BarList[elem].globalDOF[col]))
                            {
                                continue;
                            }

                            lock (this)
                            {
                                globalK.At(Codend.BarList[elem].globalDOF[row],
                                                 Codend.BarList[elem].globalDOF[col],
                                                 Codend.BarList[elem].localK[row, col]);
                            }

                        }
                    }
                });
            }

            /* ONES FOR FIXED DOF */
            Parallel.For(0, fixedDOFcount, i =>
            {
                lock (this)
                {
                    globalK.At(fixedDOF[i], fixedDOF[i], 1);
                }
            });

            /* ADDED DIAGONAL STIFFNESS*/
            Parallel.For(0, totalDOFcount, i =>
            {
                lock (this)
                {
                    globalK.At(i, i, kDiag);
                }
            });
        }

        private void UpdateShape(double[] h, double lambda)
        {
            for (int i = 0; i < totalDOFcount; i++)
            {
                if (i % 3 == 0)
                {
                    Codend.NodeList[i / 3].X += lambda * h[i];
                }

                if (i % 3 == 1)
                {
                    Codend.NodeList[i / 3].Y += lambda * h[i];
                }

                if (i % 3 == 2)
                {
                    Codend.NodeList[i / 3].Z += lambda * h[i];
                }
            }

            Codend.UpdateTotalLength();

            for (int i = 0; i < Codend.TriangleList.Count(); i++)
            {
                Codend.TriangleList[i].UpdateElement();
            }

            if (Codend.BarList.Count() != 0)
            {
                for (int i = 0; i < Codend.BarList.Count(); i++)
                {
                    Codend.BarList[i].UpdateElement();
                }
            }

        }

        private void UpdateShapeWithLS(double[] h)
        {
            double lambda = 1;                                  // new linesearch step
            double lambdac = 1;                                 // current linesearch step    
            double lambdap = 1;                                 // previous linesearch step

            int lineiter = 0;                                   // Linesearch iteration count  
            double sigma = SolverSettings.ReductionCondition;

            double[] R2hist = new double[3];
            double R0 = UpdateForceResidual();
            double R = 0;
            R2hist[0] = Math.Pow(R0, 2);     // save current residual squared to history

            for (int j = 0; j < SolverSettings.LineSearchMax; j++)
            {
                lineiter = j + 1;

                // choose lambda
                if (j == 0)                             // initially try full step
                {
                    lambda = 1;
                }
                else if (j == 1)                        // try half step if full step is not successful
                {
                    lambda = 0.5;

                }
                else if (j == SolverSettings.LineSearchMax - 1)        // If none of the method succeeded use step limit method                     
                {
                    if (h.Max() > SolverSettings.AlphaMax)
                    {
                        lambda = SolverSettings.AlphaMax / h.Max();
                    }
                    else
                    {
                        lambda = 1;
                    }
                }
                else                                    // after 2nd unsuccessful iteration try parabolic model
                {
                    lambda = ParaInterp3p(lambdac, lambdap, R2hist[2], R2hist[1], R2hist[0]);
                }

                UpdateShape(h, lambda);   // Calculate a candidate for new X with lambda step size
                AssembleForces();
                R = UpdateForceResidual();                  // Calculate residual

                if (SolverSettings.ShowLineSearchSteps)
                {
                    Console.WriteLine("{0,20} fun. eval.:{1:N0}," +
                              " Step size (lambda): {2:F3}," +
                              " Candidate residual: {3:E3}",
                              "", lineiter, lambda, R);
                }


                if (R < (1 - sigma) * R0)
                {
                    if (SolverSettings.ShowLineSearchSteps)
                    {
                        Console.WriteLine("{0,20}{1:E3} < {2:E3}:" +
                                          " Decrease condition satisfied\n",
                                          "", R, (1 - sigma) * R0);
                    }
                    break;
                }

                if (j == SolverSettings.LineSearchMax - 1)
                {
                    Console.WriteLine("{0,20}{1:E3} > {2:E3}:" +
                                      " Finishing without decrease\n",
                                     "", R, (1 - sigma) * R0);
                    break;
                }

                    UpdateShape(h, -lambda);          // set the domain geometry to initial state
                    lambdap = lambdac;                // current lambda becomes previous
                    lambdac = lambda;

                    //R2hist[0] = R2hist[1];
                    R2hist[1] = R2hist[2];            // current residual squared becomes previous
                    R2hist[2] = Math.Pow(R, 2);       // candidate residual squared becomes current
            }

        }

        private double ParaInterp3p(double lc, double lp, double fc, double fp, double f0)
        {
            /* Use line search to find out lambda that decreases the residual R

            X = X + lambda*h

            If full step (lambda = 1) and half step (lambda = 0.5) are rejected 
            due to insufficiend decrease in comparison with the residual from previous NR iteration
            next line search iterations calculate lambda using 3-point parabolic fitting with safeguard
            In addition, if pre-last iteration is neither succesful, the final step is chosen according 
            to step limit method (either full step or a fraction of a twine length depending on size of vector h)
            */

            /*      input:
            lc = current steplength
            lp = previous steplength

            fc = value of || F(x_c + lc d) ||^ 2
            fp = value of || F(x_c + lp d) ||^ 2
            f0 = value of || F(x_c)        ||^ 2
             */
            const double sigma0 = 0.1;
            const double sigma1 = 0.5;
            /*
            Compute coefficients of interpolation polynomial
            p(l) = ff0 + (c1 l + c2 l ^ 2) / d1
            p(l)' = (c1 + 2 c2 l) / d1
            p(l)'' = 2 c2 / d1

            where:
            c1 = lc*lc*(fp-f0) - lp*lp*(fc-f0);
            c2 = lp*(fc-f0) - lc*(fp-f0);
            d1 = (lc - lp) * lc * lp < 0

            so if c2 > 0 then p(l)'' < 0 we have negative curvature
            and default to lp = sigma1 * l
            */

            double c2 = lp * (fc - f0) - lc * (fp - f0);

            if (c2 >= 0)
            {
                return sigma1 * lc;
            }

            double c1 = lc * lc * (fp - f0) - lp * lp * (fc - f0);
            double lmin = -c1 * 0.5 / c2;

            /* safeguard */

            if (lmin < sigma0 * lc)
            {
                lmin = sigma0 * lc;
            }

            if (lmin > sigma1 * lc)
            {
                lmin = sigma1 * lc;
            }

            return lmin;
        }

        public void Solve()
        {
            Stopwatch Stopwatch = new Stopwatch();
            
            double[] initialX = Codend.NodeList.Select(n => n.X).ToArray();
            double[] initialY = Codend.NodeList.Select(n => n.Y).ToArray();
            double[] initialZ = Codend.NodeList.Select(n => n.Z).ToArray();

            int Iter = 0;
            int restartsCount = 0;
            bool HasConverged = false;
            bool HasDiverged = false;
            bool correctSolution = false;                       // correctness of solution flag

            int IterMax = SolverSettings.IterMax;

            double hNorm = 1;
            double R = 1;
            double lambda = 1;

            double addStiff = SolverSettings.DiagStiffness;
            double addStiff0 = SolverSettings.DiagStiffness;
            double tolStiff = SolverSettings.StiffnessTol;

            Stopwatch.Start();

            //===========================
            // added stiffness control loop
            //===========================

            while (!correctSolution)
            {
                InitializeFEM();
                double[] h = new double[totalDOFcount];

                if (restartsCount > 0)  // if the scheme doesnt work without extra stiffness  
                {
                    Iter = 0;

                    for (int i = 0; i < Codend.NodeList.Count; i++)
                    {
                        Codend.NodeList[i].X = initialX[i];
                        Codend.NodeList[i].Y = initialY[i];
                        Codend.NodeList[i].Z = initialZ[i];
                    }
                    
                    if (addStiff0 == 0)
                    {
                        addStiff0 = 1;
                    }
                }

                if (restartsCount > 100)
                {
                    Console.WriteLine("Stopped, due to too many restarts");
                    break;
                }
                addStiff = addStiff0;
                tolStiff = SolverSettings.StiffnessTol;

                UpdateShape(h, 0);
                AssembleForcesAndStiffness(addStiff); // here only forces actually

                R = UpdateForceResidual();

                PrintIter(Iter, 0, addStiff, R, hNorm);

                //===========================
                // Newton-Raphson loop
                //===========================

                while (!HasConverged)
                {
                    // check if diverged
                    if (R > SolverSettings.ResidualMax)
                    {
                        HasDiverged = true;
                        break;
                    }

                    if (Iter == IterMax)
                    {
                        HasDiverged = true;
                        break;
                    }
                    else
                    {
                        Iter++;
                    }

                    if (R < tolStiff)
                    {
                        tolStiff *= SolverSettings.ReduceStiffnessBy;
                        addStiff *= SolverSettings.ReduceStiffnessBy;
                    }

                    if (R < SolverSettings.ResidualTol) //hNorm && < SolverSettings.DisplacementTol)
                    {
                        HasConverged = true;
                    }

                    h = SolveSparseSystem(h);

                    if (SolverSettings.IncludeLineSearch)
                    {
                        UpdateShapeWithLS(h);
                    }
                    else
                    {
                        lambda = InexactLineSearch();
                        UpdateShape(h, lambda);
                    }
                   
                    hNorm = DisplacmentNorm(h);

                    Codend.UpdateTotalLength();
                    AssembleForcesAndStiffness(addStiff);
                    R = UpdateForceResidual();

                    if (HasConverged)
                    {
                        Console.WriteLine("\nConvergence acheived in {0} iterations and {1} restarts in {2:mm\\:ss} [min:sec]\n",
                                Iter, restartsCount, Stopwatch.Elapsed);
                    }
                    else
                    {
                        PrintIter(Iter, 0, addStiff, R, hNorm);
                    }

                }

                //===================
                // Restart routine               
                //===================

                if (HasDiverged)
                {
                    restartsCount += 1;
                    addStiff0 *= SolverSettings.IncreaseStiffnessBy;
                    HasDiverged = false;
                    Console.WriteLine("\nIteration diverges or the algebraic system is too ill-conditioned." +
                                      "\nRestarting with higher added diagonal stiffness.\n");
                }
                else
                {
                    correctSolution = true; // correct solution found
                }                
            }  
            Stopwatch.Stop();
        }

        private double InexactLineSearch()
        {
            return 1;
        }

        private void PrintIter(int iter, int funEval, double addStiff, double R, double hNorm)
        {
            if (SolverSettings.IncludeLineSearch)
            {
                Console.WriteLine(String.Format("Iter: {0,-7:D}" +
                                            "LS steps: {1,-6:D}" +
                                            "addStiff: {2,-10:e2}{3,-8:S}" +
                                            "R: {4,-10:e2}{5,-6:S}" +
                                            "|h|: {6,-10:e2}{7,-6:S}", iter, funEval, addStiff, "N/m", R, "N", hNorm, "m"));
            }
            else
            {
                Console.WriteLine(String.Format("Iter: {0,-7:D}" +
                            "addStiff: {1,-10:e2}{2,-8:S}" +
                            "R: {3,-10:e2}{4,-6:S}" +
                            "|h|: {5,-10:e2}{6,-6:S}", iter, addStiff, "N/m", R, "N", hNorm, "m"));
            }

        }

        private double[] SolveSparseSystem(double[] h)
        {
            var order = ColumnOrdering.MinimumDegreeAtPlusA;
            double tolerance = 1.0;
            var lu = SparseLU.Create(Converter.ToCompressedColumnStorage(globalK), order, tolerance);
            lu.Solve(globalF, h);
            return h;
        }

        public void Simulate()
        {
            //if (!File.Exists(FilePaths.Output3dResults))
            //{
                StartResultsFile(FilePaths.Output3dResults);
            //}
           
            for (int i = 0; i < Catch.Count; i++)
            {
                if (SolverSettings.InitialShape.Contains("Axi"))// && Codend.SameLengthPanels())
                {
                    Codend.SetAxiInitialShape(FilePaths, Catch, i, Towing);
                }
                else
                {
                    Codend.SetCylInitialShape();
                }

                Codend.RestrainEntrance();
                //CalculateBoundaryDOF();

                ApplyTowing();
                Catch.Apply(Codend, i);

                Solve();

                //SaveFemShape(FilePaths.Output3dShapes, (int)Catch.BlockedMeshes[i]);
                SaveFemShape(FilePaths.Output3dShapes, i);
                AppendResults(FilePaths.Output3dResults);
            }
        }

        public void ConvergenceStudy()
        {
            StartResultsFile(FilePaths.Output3dResults);
            string dofCount = "";

            for (int i = 40; i >= 0; i--)
            {
                Codend.MeshSettings.ElemAcrossPanel = 10 + i;
                Codend.MeshSettings.ElemAlongPanel = 10 + i;

                if (SolverSettings.InitialShape.Contains("Axi") && Codend.SameLengthPanels())
                {
                    Codend.SetAxiInitialShape(FilePaths, Catch, i, Towing);
                }
                else
                {
                    Codend.SetCylInitialShape();
                }

                Codend.RestrainEntrance();
                CalculateBoundaryDOF();

                ApplyTowing();
                Catch.Apply(Codend, 0);

                Solve();
                dofCount += String.Format(" {0}",totalDOFcount);

                AppendResults(FilePaths.Output3dResults);
            }
            Console.WriteLine(dofCount);
        }

        /* results */

        public double MaxRadius()
        {
            double rMax = Codend.EntranceRadius;                           // maximum radius
            double rTest = 0;
            for (int i = 0; i < Codend.NodeList.Count; i++)
            {
                rTest = Math.Sqrt(Math.Pow(Codend.NodeList[i].X, 2) + Math.Pow(Codend.NodeList[i].Z, 2));
                if (rMax < rTest)
                {
                    rMax = rTest;
                }
            }
            return rMax;
        }

        public double EntranceDrag()
        {
            double ReY = 0; // Take Y reaction

            Parallel.For(0, Codend.TriangleList.Count, elem =>
            {
                Codend.TriangleList[elem].UpdateElementTotalForce();

                for (int row = 0; row < 9; row++)
                {
                    if (!IsFreeDOF(Codend.TriangleList[elem].globalDOF[row]) && row % 3 == 1)
                    {
                        lock (this)
                        {
                            ReY += Codend.TriangleList[elem].localF[row];
                        }
                    }
                }
            });

            return ReY;   // resultant entrance reaction
        }

        public List<Node> CatchBound()
        {
            double yMin = Codend.Length;
            foreach (var t in Codend.TriangleList)
            {
                if(t.HasCatch)
                {
                    yMin = (yMin > t.n1.Y) ? t.n1.Y : yMin;
                    yMin = (yMin > t.n2.Y) ? t.n2.Y : yMin;
                    yMin = (yMin > t.n3.Y) ? t.n3.Y : yMin;
                }
            }

            return Codend.NodeList.FindAll(n => Math.Abs(n.Y - yMin) < 1e-3);
        }
        
        public double CatchDrag()
        {
            //double rCatch = X[2 * ncp - 1];      // radius at the beginning of the catch
            return 0;//pi * Math.Pow(rCatch, 2) * P; // resultant catch drag
        }

        public double CatchThickness()
        {
            var CatchBoundNodes = CatchBound();
            double catchY = 0;
            foreach (var n in CatchBoundNodes)
            {
                catchY += n.Y;
            }
            catchY /= CatchBoundNodes.Count;

            return Codend.Length - catchY;
        }

        public double CatchVolume()
        {
            var CatchBoundNodes = CatchBound();
            var catchLim = new Node(-1,0,0,0);
            foreach (var n in CatchBoundNodes)
            {
                catchLim.X += n.X;
                catchLim.Y += n.Y;
                catchLim.Z += n.Z;
            }
            catchLim.X /= CatchBoundNodes.Count;
            catchLim.Y /= CatchBoundNodes.Count;
            catchLim.Z /= CatchBoundNodes.Count;

            double V = 0;

            for (int i = 0; i < CatchBoundNodes.Count; i++)
            {
                if(i == CatchBoundNodes.Count - 1)
                {
                    V += SignedVolume(catchLim, CatchBoundNodes[0], CatchBoundNodes[i]);
                }
                else
                {
                    V += SignedVolume(catchLim, CatchBoundNodes[i + 1], CatchBoundNodes[i]);
                }
            }
            foreach (var t in Codend.TriangleList)
            {
                if (t.HasCatch)
                {
                    V += SignedVolume(t.n1, t.n2, t.n3);
                }

            }
            return V;
        }

        public double CatchSurface()
        {
            double S = 0;
            foreach (var t in Codend.TriangleList)
            {
                if (t.HasCatch)
                {
                    S += 0.5 * Math.Abs(t.n1.X * (t.n2.Y - t.n3.Y) + 
                                        t.n2.X * (t.n3.Y - t.n1.Y) + 
                                        t.n3.X * (t.n1.Y - t.n2.Y));
                }
            }
            return S;
        }

        private double SignedVolume(Node n1, Node n2, Node n3)
        {
            return 
                (n1.X * n2.Y * n3.Z +
                -n1.X * n3.Y * n2.Z +
                -n2.X * n1.Y * n3.Z +
                 n2.X * n3.Y * n1.Z +
                 n3.X * n1.Y * n2.Z +
                -n3.X * n2.Y * n1.Z) / 6.0;
        }

        private double TestSignedVolume()
        {
            
            var pts = new List<Node>()
            {
                new Node(0,-0.5,1,-0.5),
                new Node(1, 0.5,1,-0.5),
                new Node(2, 0.5,1, 0.5),
                new Node(3,-0.5,1, 0.5),
                new Node(4,-5,3,-5),
                new Node(5, 5,3,-5),
                new Node(6, 5,3, 5),
                new Node(7,-5,3, 5)
            };

            int[,] T = new int[,] {
                {0, 1, 2 },
                {0, 2, 3 },
                {1, 5, 6 },
                {1, 6, 2 },
                {0, 3, 7 },
                {0, 7, 4 },
                {2, 6, 7 },
                {2, 7, 3 },
                {1, 0, 4 },
                {1, 4, 5 },
                {5, 4, 7 },
                {5, 7, 6 } };

            double V = 0;
            for (int i = 0; i < T.GetLength(0); i++)
            {
                V += SignedVolume(pts[T[i,0]], pts[T[i, 1]], pts[T[i, 2]]);
            }
            return V;
        }

        /* output */

        private void StartResultsFile(string path)
        {
            using (StreamWriter ResultFile = new StreamWriter(path))
            {
                ResultFile.WriteLine("{0,-20}{1,-20}{2,-20}{3,-20}{4,-20}{5,-20}",
                                     "Length",
                                     "Max radius",
                                     "Catch thickness",
                                     "Catch surface",
                                     "Catch volume",
                                     "Total reaction");
            }
        }

        private void AppendResults(string path)
        {
            using (StreamWriter ResultFile = File.AppendText(path))
            {
                ResultFile.WriteLine("{0,-20:E5}{1,-20:E5}{2,-20:E5}{3,-20:E5}{4,-20:E5}{5,-20:E5}",
                     Codend.Length,
                     MaxRadius(),
                     CatchThickness(),
                     CatchSurface(),
                     CatchVolume(),
                     EntranceDrag());
            }
        }

        private void SaveFemShape(string path, int label)
        {
            string extension = ".txt";           
            string pathWithLabel = path.Substring(0, path.Length - extension.Length) + label + extension;
            using (StreamWriter ResultFile = new StreamWriter(pathWithLabel))
            {
                /*NODES*/
                ResultFile.WriteLine("NODES {0}", Codend.NodeList.Count);
                ResultFile.WriteLine("{0,-5}{1,-15}{2,-15}{3,-15}",
                                     "ID", "X [m]", "Y [m]", "Z [m]");

                foreach (var node in Codend.NodeList)
                {
                    ResultFile.WriteLine("{0,-5:D}{1,-15:E5}{2,-15:E5}{3,-15:E5}", node.ID, node.X, node.Y, node.Z);
                }

                /*TRIANGLES*/
                ResultFile.WriteLine("TRIANGLES {0}", Codend.TriangleList.Count);
                ResultFile.WriteLine("{0,-10}{1,-10}{2,-10}{3,-10}{4,-10}{5,-10}{6,-10}{7,-10}{8,-10}{9,-10}{10,-10}",
                                      "ID", "n1 ID", "n2 ID", "n3 ID", "n1 u", "n1 v", "n2 u", "n2 v", "n3 u", "n3 v", "Has Catch");

                foreach (var tri in Codend.TriangleList)
                {
                    ResultFile.WriteLine("{0,-10:D}{1,-10:D}{2,-10:D}{3,-10:D}{4,-10:F2}{5,-10:F2}{6,-10:F2}{7,-10:F2}{8,-10:F2}{9,-10:F2}{10,-10:D}",
                                       tri.ID, tri.n1.ID, tri.n2.ID, tri.n3.ID,
                                       tri.uCoord[0], tri.vCoord[0], tri.uCoord[1], tri.vCoord[1], tri.uCoord[2], tri.vCoord[2],
                                       tri.HasCatch);
                }

                /*BARS*/
                ResultFile.WriteLine("BARS {0}", Codend.BarList.Count);
                ResultFile.WriteLine("{0,-10}{1,-10}{2,-10}{3,-10}",
                                     "ID", "n1 ID", "n2 ID", "L [m]");

                foreach (var bar in Codend.BarList)
                {
                    ResultFile.WriteLine("{0,-10:D}{1,-10:D}{2,-10:D}{3,-15:E5}",
                                       bar.ID, bar.n1.ID, bar.n2.ID, bar.L);

                }
            }
        }
    }
}

