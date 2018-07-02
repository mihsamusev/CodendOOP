using System;
using System.Collections.Generic;
using System.Reflection;
using System.IO;
using System.Globalization;

namespace CodendOOP
{
    class Program
    {
        static void Main(string[] args)
        {
            CultureInfo.DefaultThreadCurrentCulture = new CultureInfo("en-US");
            try
            {
                Console.WindowWidth = 100;
                Console.WindowHeight = 50;
            }
            catch
            {

            }
            Console.WindowHeight = 50;

            //=========================
            // INPUT
            //=========================

            var inputReader = new InputReader();
            inputReader.PrintInfo();

            var nPanels = inputReader.CountPanels();
            var nPanelMaterials = inputReader.CountPanelMaterials();

            //=========================================
            // INITIALIZE PANEL MATERIALS
            //=========================================

            List<PanelMaterial> PanelMaterialList = new List<PanelMaterial>();
            for (int i = 0; i < nPanelMaterials; i++)
            {
                PanelMaterialList.Add(new PanelMaterial(i + 1, inputReader));
            }

            //foreach (var m in PanelMaterialList)
            //{
            //    m.PrintInfo();
            //}

            //=========================================
            // INITIALIZE PANELS
            //=========================================

            List<Panel> PanelList = new List<Panel>();

            for (int i = 0; i < nPanels; i++)
            {
                var panelInput = inputReader.ReadPanelInput(i + 1);

                if (panelInput.type.Equals("DiamondPanel", StringComparison.InvariantCultureIgnoreCase))
                {
                    PanelList.Add( new DiamondMeshPanel(
                        LengthInMeshes: panelInput.meshesAlong,
                        WidthInMeshes: panelInput.meshesAcross,
                        Orientation: panelInput.orientation,
                        Material: PanelMaterialList[panelInput.materialID - 1]));
                }
                else if (panelInput.type.Equals("SquarePanel", StringComparison.InvariantCultureIgnoreCase))
                {
                    PanelList.Add(new SquareMeshPanel(
                        LengthInMeshes: panelInput.meshesAlong,
                        WidthInMeshes: panelInput.meshesAcross,
                        Orientation: panelInput.orientation,
                        Material: PanelMaterialList[panelInput.materialID - 1]));
                }
                else
                {
                    throw new ArgumentException("Panel types should be ither \'DiamondPanel\' or \'SquarePanel\'");
                }
            }

            //foreach (var p in PanelList)
            //{
            //    p.PrintInfo();
            //}

            //=========================================
            // INITIALIZE SELVEDGES
            //=========================================

            List<Selvedge> SelvedgeList = new List<Selvedge>();

            if (inputReader.SelvedgesIncluded())
            {
                for (int i = 0; i < nPanels; i++)
                {
                    SelvedgeList.Add(new Selvedge(i, inputReader));
                }
            }

            //=========================================
            // INITIALIZE ROUND STRAPS
            //=========================================

            List<RoundStrap> RoundStrapList = new List<RoundStrap>();

            if (inputReader.RoundStrapsIncluded())
            {
                int strapCount = inputReader.CountRoundStraps();
                for (int i = 0; i < strapCount; i++)
                {
                    RoundStrapList.Add(new RoundStrap(i + 1, inputReader));
                }
            }

            //=========================================
            // INITIALIZE CODEND
            //=========================================

            double R = inputReader.ReadCodendEntranceRadius();
            var Codend = new Codend(PanelList, R)
            {
                SelvedgeList = SelvedgeList,
                StrapList = RoundStrapList,
                MeshSettings = new MeshSettings(inputReader)
            };

            //=========================================
            // INITIALIZE CATCH
            //=========================================

            var Catch = new Catch(inputReader);

            //=========================================
            // INITIALIZE TOWING
            //=========================================

            var TowingAlongY = new Towing(inputReader);

            //=========================================
            // INITIALIZE SIMULATION
            //=========================================

            var TowingSimulation = new Simulation(Codend, Catch, TowingAlongY)
            {
                SolverSettings = new SolverSettings(inputReader),
                FilePaths = inputReader
            };


            TowingSimulation.Simulate();

            ////TowingSimulation.ConvergenceStudy();

            //Console.ReadKey();
        }

        //======================================================================
        // HELPER METHODS
        //======================================================================

        static string AssemblyPath()
        {
            string fullPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().GetName().CodeBase);
            return fullPath.Replace("file:\\", "");
        }

        static void PrintVector(double[] Vec)
        {
            for (int i = 0; i < Vec.Length; i++)
            {
                Console.WriteLine("{0,10:F4}", Vec[i]);
            }
        }

        static void PrintMatrix(double[,] Mat)
        {
            for (int i = 0; i < Mat.GetLength(0); i++)
            {
                for (int j = 0; j < Mat.GetLength(1); j++)
                {
                    Console.Write("{0,15:F3}", Mat[i, j]);
                }
                Console.Write("\n");
            }

        }

        static void SavePanel(List<Node> NodeList, List<TriangleElement> TriangleList, string path)
        {
            using (StreamWriter ResultFile = new StreamWriter(path))
            {
                /*NODES*/
                ResultFile.WriteLine("NODES {0}", NodeList.Count);
                ResultFile.WriteLine("{0,-5}{1,-15}{2,-15}{3,-15}",
                                     "ID", "X [m]", "Y [m]", "Z [m]");

                foreach (var node in NodeList)
                {
                    ResultFile.WriteLine("{0,-5:D}{1,-15:E5}{2,-15:E5}{3,-15:E5}", node.ID, node.X, node.Y, node.Z);
                }

                /*TRIANGLES*/
                ResultFile.WriteLine("TRIANGLES {0}", TriangleList.Count);
                ResultFile.WriteLine("{0,-10}{1,-10}{2,-10}{3,-10}{4,-10}{5,-10}{6,-10}{7,-10}{8,-10}{9,-10}{10,-10}",
                                      "ID", "n1 ID", "n2 ID", "n3 ID", "n1 u", "n1 v", "n2 u", "n2 v", "n3 u", "n3 v", "Has Catch");

                foreach (var tri in TriangleList)
                {
                    ResultFile.WriteLine("{0,-10:D}{1,-10:D}{2,-10:D}{3,-10:D}{4,-10:F2}{5,-10:F2}{6,-10:F2}{7,-10:F2}{8,-10:F2}{9,-10:F2}{10,-10:D}",
                                       tri.ID, tri.n1.ID, tri.n2.ID, tri.n3.ID,
                                       tri.uCoord[0], tri.vCoord[0], tri.uCoord[1], tri.vCoord[1], tri.uCoord[2], tri.vCoord[2],
                                       tri.HasCatch);
                }
            }
        }

    }
}

