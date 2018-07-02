using System;
using System.IO;
using System.Linq;

namespace CodendOOP
{
    class Catch
    {
        //=========================
        // variables
        //=========================

        public double[] BlockingRatio;
        public double[] BlockedMeshes;
        public string applyMethod;
        public int Count;
        private readonly string method1 = "ByBlockedMeshes";
        private readonly string method2 = "ByBlockingRatio";
        //public double DragCoefficient = 1.4;

        //=========================
        // constructors
        //=========================

        public Catch(double[] BlockedMeshes)
        {
            applyMethod = "ByBlockedMeshes";
            this.BlockedMeshes = BlockedMeshes;
            Count = BlockedMeshes.Length;
        }

        public Catch(InputReader filePath)
        {
            LoadInput(filePath.Input3d);
        }

        //=========================
        // methods
        //=========================

        public void ApplyByBlockingRatio(TriangleElement Tri, double CodendLength, double LengthRatio)
        {
            double TriMinY = (new double[3] { Tri.n1.Y, Tri.n2.Y, Tri.n3.Y }).Min();
            if (TriMinY >= (1 - LengthRatio) * CodendLength)
            {
                Tri.HasCatch = true;
            }
            else
            {
                Tri.HasCatch = false;
            }
        }

        public void ApplyToDiamondByBlockedMeshes(TriangleElement Tri, double CodendLengthInMeshes, double NumBlockedMeshes)
        {
            double vTriMin = (Tri.vCoord).Min();
            if (vTriMin >= (CodendLengthInMeshes - NumBlockedMeshes))
            {
                Tri.HasCatch = true;
            }
            else
            {
                Tri.HasCatch = false;
            }
        }

        public void ApplyToSquareByBlockedMeshes(TriangleElement Tri, double CodendLengthInMeshes, double NumBlockedMeshes)
        {
            double VTriMin = (Tri.VCoord).Min();
            if (VTriMin >= (CodendLengthInMeshes - NumBlockedMeshes))
            {
                Tri.HasCatch = true;
            }
            else
            {
                Tri.HasCatch = false;
            }
        }

        //public void LabelCatchBoundary(NetTriangle Tri, int label)
        //{
        //    if (Tri.HasCatch)
        //    {
        //        double min = new double[3] { Tri.n1.Y, Tri.n2.Y, Tri.n3.Y }.Min();
        //        if(Tri.n1.Y == min)
        //        {
        //            Tri.n1.Label = label;
        //        }
        //        if (Tri.n2.Y == min)
        //        {
        //            Tri.n2.Label = label;
        //        }
        //        if (Tri.n3.Y == min)
        //        {
        //            Tri.n3.Label = label;
        //        }
        //    }
        //}

        public void Apply(Codend Codend, int catchNumber)
        {
            int nTri = Codend.TriangleList.Count;

            if(nTri == 0)
            {
                throw new ArgumentException("Codend should be meshed in order to apply catch. Provide initial shape");
            }

            // Apply by blocked meshes
            if (applyMethod.Equals(method1, StringComparison.InvariantCultureIgnoreCase))
            {
                if (Codend.SameLengthPanels())
                {
                    double CodendLengthInMeshes = Codend.PanelList[0].LengthInMeshes;
                    double NumBlockedMeshes = BlockedMeshes[catchNumber];


                    if (Codend.PanelList[0].Material.MeshType.Equals("Diamond", StringComparison.InvariantCultureIgnoreCase))
                    {
                        for (int i = 0; i < nTri; i++)
                        {
                            ApplyToDiamondByBlockedMeshes(Codend.TriangleList[i], CodendLengthInMeshes, NumBlockedMeshes);
                        }
                    }
                    else if (Codend.PanelList[0].Material.MeshType.Equals("Square", StringComparison.InvariantCultureIgnoreCase))
                    {
                        for (int i = 0; i < nTri; i++)
                        {
                            ApplyToSquareByBlockedMeshes(Codend.TriangleList[i], CodendLengthInMeshes, NumBlockedMeshes);
                        }
                    }
                }
                else
                {
                    throw new ArgumentException("Panels must have same amount of meshes along to apply catch by blocked meshes method");
                }
            }
            // Apply by length ratios
            else if (applyMethod.Equals(method2, StringComparison.InvariantCultureIgnoreCase))
            {
                if (BlockingRatio[catchNumber] != 0)
                {
                    double blockingRatio = BlockingRatio[catchNumber];
                    for (int i = 0; i < nTri; i++)
                    {
                        ApplyByBlockingRatio(Codend.TriangleList[i], Codend.Length, blockingRatio);
                    }
                }
                else
                {
                    throw new ArgumentException("Blocking ration cannot be 0 to apply catch by blocking ratio method");
                }
            }           
        }

        private void LoadInput(string inputPath)
        {
            string targetWord = "CatchCount";
            string[] lines = File.ReadAllLines(inputPath);
            int currentLine;
            string[] parts;
            bool catchFound = false;

            currentLine = 0;
            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    catchFound = true;
                    break;
                }
                currentLine++;
            }

            if (!catchFound)
            {
                throw new IOException("Wrong input file format, catch is not found!");
            }
            else
            {
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                Count = Convert.ToInt32(parts[1]);

                currentLine++;
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                if (parts[0].Equals("ApplyMethod", StringComparison.InvariantCultureIgnoreCase))
                {
                    applyMethod = parts[1];
                }
                else
                {
                    throw new IOException("\'ApplyMethod\' is not defined for catch input!");
                }



                if (Count > 0 && applyMethod.Equals(method1, StringComparison.InvariantCultureIgnoreCase))
                {
                    currentLine++;
                    BlockedMeshes = new double[Count];
                    for (int i = 0; i < Count; i++)
                    {
                        parts = lines[currentLine + i].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        BlockedMeshes[i] = Convert.ToDouble(parts[0]);
                    }
                }
                else if (Count > 0 && applyMethod.Equals(method2, StringComparison.InvariantCultureIgnoreCase))
                {
                    double ratio;
                    currentLine++;
                    BlockingRatio = new double[Count];
                    for (int i = 0; i < Count; i++)
                    {
                        parts = lines[currentLine + i].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        ratio = Convert.ToDouble(parts[0]);
                        if (ratio > 0 && ratio <= 1)
                        {
                            BlockingRatio[i] = ratio;
                        }
                        else
                        {
                            throw new IOException("Blocking ratios should be between 0 and 1");
                        }                                                  
                    }
                }
            }
        }



    }
}
