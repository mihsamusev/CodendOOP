using System;
using System.Collections.Generic;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace CodendOOP
{
    class Codend
    {
        //==================
        // fields
        //==================

        /*real life object lists*/
        public List<Panel> PanelList;
        public List<Selvedge> SelvedgeList { get; set; }
        public List<RoundStrap> StrapList { get; set; }

        /*finite element object lists*/
        public List<Node> NodeList = new List<Node>();
        public List<TriangleElement> TriangleList = new List<TriangleElement>();
        public List<Restraint> RestraintList = new List<Restraint>();
        public List<BarElement> BarList = new List<BarElement>();

        public readonly int PanelCount;
        public double Length;
        public double EntranceRadius;

        public MeshSettings MeshSettings;

        //==================
        // constructor
        //==================

        public Codend(List<Panel> PanelList, double EntranceRadius)
        {
            this.PanelList = PanelList;
            PanelCount = PanelList.Count;
            this.EntranceRadius = EntranceRadius;
            MeshSettings = new MeshSettings();
        }

        //==================
        // methods
        //==================

        public void SetCylInitialShape()
        {
            InitializeFemLists();

            // find maximum panel width and height
            double maxLength = PanelList[0].Length;
            double maxWidth = PanelList[0].Width;

            for (int i = 1; i < PanelCount; i++)
            {
                if (maxLength < PanelList[i].Length)
                {
                    maxLength = PanelList[i].Length;
                }
                if (maxWidth < PanelList[i].Width)
                {
                    maxWidth = PanelList[i].Width;
                }
            }

            double totalAngle = 2 * Math.PI / PanelCount;

            for (int i = 0; i < PanelCount; i++)
            {
                PanelList[i].CreateUnitMesh(MeshSettings);

                if (StrapList.Count != 0)
                {
                    for (int j = 0; j < StrapList.Count; j++)
                    {
                        PanelList[i].LabelNodesForStrap(StrapList[j]);
                    }
                }

                PanelList[i].ScaleUnitMesh(maxLength, maxWidth);
                PanelList[i].MapToCyllinder(EntranceRadius, -totalAngle, totalAngle / 2 - i * totalAngle);
            }

            if (SelvedgeList.Count != 0)
            {
                CreateAllSelvedges();
            }

            if (StrapList.Count != 0)
            {
                CreateAllStraps();
            }

            JoinMeshedPanels();

            CloseCylCodend(); // set last and pre last x and z to 0.5 and 0 of their values
          
            UpdateTotalLength();
        }

        private void CreateAllSelvedges()
        {
            for (int i = 0; i < PanelCount; i++)
            {
                SelvedgeList[i].SetNodeList(PanelList[i].LeftSelvedgeNodes);
                SelvedgeList[i].CreateBarElements();
            }
        }

        private void CreateAllStraps()
        {
            List<Node> strapNodes;
            for (int i = 0; i < StrapList.Count; i++)
            {
                strapNodes = new List<Node>();
                for (int j = 0; j < PanelCount; j++)
                {
                    strapNodes.AddRange(PanelList[j].NodeList.FindAll(n => n.Label == StrapList[i].Label));
                    
                }
                StrapList[i].SetNodeList(strapNodes);
                StrapList[i].CreateBarElements();
            }
        }

        private void CloseCylCodend()
        {
            for (int i = 0; i < NodeList.Count; i++)
            {
                if (NodeList[i].Label == 2) // squeeze pre-last ring of points
                {
                    NodeList[i].X = 0.5 * NodeList[i].X;
                    NodeList[i].Z = 0.5 * NodeList[i].Z;
                }

                if (NodeList[i].Label == 1) // squeeze last ring of points
                {
                    NodeList[i].X = 0;
                    NodeList[i].Z = 0;
                }
            }
            UpdateGlobalTopology();
        }

        public void SetAxiInitialShape(InputReader filePaths, Catch Catch, int catchNumber, Towing towing)
        {
            InitializeFemLists();

            if (Catch.applyMethod.Equals("ByBlockedMeshes", StringComparison.InvariantCultureIgnoreCase))
            {
                WriteAxiInput(filePaths, PanelList[0].Material, Catch.BlockedMeshes[catchNumber], towing.Speed);
            }
            else if (Catch.applyMethod.Equals("ByBlockingRatio", StringComparison.InvariantCultureIgnoreCase))
            {
                double ratio = Catch.BlockingRatio[catchNumber];
                double blockedMeshes = Math.Floor(PanelList[0].LengthInMeshes * ratio);
                WriteAxiInput(filePaths, PanelList[0].Material, blockedMeshes, towing.Speed);
            }
            else
            {
                throw new ArgumentException("Unknown catch application method");
            }

            Process p = new Process();
            //p.StartInfo.WindowStyle = ProcessWindowStyle.Hidden;
            p.StartInfo.FileName = filePaths.AxiExe;
            p.Start();
            p.WaitForExit();

            double[] AxiX = RemoveAxiKnots(ReadAxiOutput(filePaths.OutputAxiShapes));

            double totalAngle = 2 * Math.PI / PanelCount;

            for (int i = 0; i < PanelCount; i++)
            {
                PanelList[i].CreateUnitMesh(MeshSettings);

                if (StrapList.Count != 0)
                {
                    for (int j = 0; j < StrapList.Count; j++)
                    {
                        PanelList[i].LabelNodesForStrap(StrapList[j]);
                    }
                }

                PanelList[i].MapToAxiModel(AxiX, -totalAngle, totalAngle / 2 - i * totalAngle);
            }

            if (SelvedgeList.Count != 0)
            {
                CreateAllSelvedges();
            }

            if (StrapList.Count != 0)
            {
                CreateAllStraps();
            }

            JoinMeshedPanels();

            UpdateTotalLength();
        }

        public void UpdateGlobalTopology()
        {
            int count = 0;
            for (int i = 0; i < NodeList.Count; i++)
            {
                if (NodeList[i].CopyOf != -1)
                {
                    NodeList[i].ID = NodeList[i].CopyOf;
                }
                else
                {
                    NodeList[i].ID = count;
                    count++;
                }

                for (int j = i + 1; j < NodeList.Count; j++)
                {
                    if (NodeList[i].IsEqual(NodeList[j]))
                    {
                        NodeList[j].CopyOf = NodeList[i].ID;
                    }
                }
            }

            count = 0;
            for (int i = 0; i < TriangleList.Count; i++)
            {
                TriangleList[i].UpdateNodeID(NodeList);

                if (TriangleList[i].n1 == TriangleList[i].n2 || // label singular triangles 
                    TriangleList[i].n2 == TriangleList[i].n3 ||
                    TriangleList[i].n3 == TriangleList[i].n1)
                {
                    TriangleList[i].Label = 1;
                }
                else
                {
                    TriangleList[i].ID = count;
                    count++;
                }
            }

            count = 0;
            for (int i = 0; i < BarList.Count; i++)
            {
                BarList[i].UpdateNodeID(NodeList);

                if (BarList[i].n1 == BarList[i].n2) // label singular bars

                {
                    BarList[i].Label = 1;
                }
                else
                {
                    BarList[i].ID = count;
                    count++;
                }
            }

            NodeList.RemoveAll(n => n.CopyOf != -1);    // remove all non unique nodes
            TriangleList.RemoveAll(t => t.Label == 1);  // remove all singular triangles
            BarList.RemoveAll(b => b.Label == 1);       // remove all singular bars
        }

        public void RestrainEntrance()
        {
            for (int i = 0; i < NodeList.Count; i++)
            {
                if (NodeList[i].Label == 0)
                {
                    RestraintList.Add(new Restraint(nodeID: NodeList[i].ID, xIsFixed: true, yIsFixed: true, zIsFixed: true));
                }
            }
        }

        public void RestrainEnd()
        {
            for (int i = 0; i < NodeList.Count; i++)
            {
                if (NodeList[i].Label == 1)
                {
                    RestraintList.Add(new Restraint(nodeID: NodeList[i].ID, xIsFixed: false, yIsFixed: true, zIsFixed: true));
                }
            }
        }

        public void UpdateTotalLength()
        {
            for (int i = 0; i < NodeList.Count; i++)
            {
                if (NodeList[i].Label == 1)
                {
                    Length = NodeList[i].Y;
                    break;
                }
            }

        }

        public bool SameLengthPanels()
        {
            double Meshes0 = PanelList[0].LengthInMeshes;

            for (int i = 1; i < PanelCount; i++)
            {
                if(Meshes0 != PanelList[i].LengthInMeshes)
                {
                    return false;
                }
            }

            return true;
        }



        private void InitializeFemLists()
        {
            NodeList = new List<Node>();
            TriangleList = new List<TriangleElement>();
            RestraintList = new List<Restraint>();
            BarList = new List<BarElement>();
        }

        private void JoinMeshedPanels()
        {
            for (int i = 0; i < PanelList.Count; i++)
            {
                NodeList.AddRange(PanelList[i].NodeList);
                TriangleList.AddRange(PanelList[i].TriangleList);
            }

            if (SelvedgeList.Count != 0)
            {
                for (int i = 0; i < SelvedgeList.Count; i++)
                {
                    BarList.AddRange(SelvedgeList[i].BarList);
                }
            }

            if (StrapList.Count != 0)
            {
                for (int i = 0; i < StrapList.Count; i++)
                {
                    BarList.AddRange(StrapList[i].BarList);
                }
            }

            UpdateGlobalTopology();
        }

        private void WriteAxiInput(InputReader filePaths, PanelMaterial Material, double BlockedMeshes, double TowingSpeed)
        {
            using (StreamWriter ResultFile = new StreamWriter(filePaths.InputAxi))
            {
                ResultFile.WriteLine("Paths");
                ResultFile.WriteLine("OutputShapes\t{0}", filePaths.OutputAxiShapes);
                ResultFile.WriteLine("Material");
                ResultFile.WriteLine("MeshSide\t{0}", Material.MeshSide);
                ResultFile.WriteLine("KnotSize\t{0}", Material.MeshSide / 1000);
                ResultFile.WriteLine("TwineEA\t{0}", Material.EA);
                ResultFile.WriteLine("KnotEA\t{0}", Material.EA);
                ResultFile.WriteLine("MeshOrientation\t0");
                ResultFile.WriteLine("\nCodend");
                ResultFile.WriteLine("MeshesAlong\t{0}", PanelList[0].LengthInMeshes);
                ResultFile.WriteLine("MeshesAround\t{0}", PanelCount * PanelList[0].WidthInMeshes);
                ResultFile.WriteLine("EntranceRadius\t{0}", EntranceRadius);
                ResultFile.WriteLine("\nCatch");
                ResultFile.WriteLine("Count\t1");
                ResultFile.WriteLine((int)BlockedMeshes);
                ResultFile.WriteLine("\nTowing");
                ResultFile.WriteLine("TowingSpeed\t{0}", TowingSpeed);
                ResultFile.WriteLine("DiagStiffness\t{0}", 100);
            }
        }

        private double[] ReadAxiOutput(string path)
        {
 
            string[] lines = File.ReadAllLines(path);
            double[] X = new double[lines.Length];

            for (int i = 0; i < lines.Length; i++)
            {
                X[i] = Convert.ToDouble(lines[i]);
            }

            return X;
        }

        private double[] RemoveAxiKnots(double[] X)
        {
            int nodesOld = X.Length / 2;
            int nodesNew = 1 + (nodesOld - 1) / 2;

            double[] Xnoknot = new double[2 * nodesNew];
            Xnoknot[0] = X[0];                                // very first node is unchanged
            Xnoknot[1] = X[1];
            Xnoknot[2 * nodesNew - 2] = X[2 * nodesOld - 2];  // very last node is unchanged
            Xnoknot[2 * nodesNew - 1] = X[2 * nodesOld - 1];

            for (int i = 1; i < nodesNew - 1; i++)              // mean of left and right knot node
            {
                Xnoknot[2 * i]     = 0.5 * (X[4 * i] +     X[4 * i + 2]);
                Xnoknot[2 * i + 1] = 0.5 * (X[4 * i + 1] + X[4 * i + 3]);
            }
            return Xnoknot;
        }

    }
}
