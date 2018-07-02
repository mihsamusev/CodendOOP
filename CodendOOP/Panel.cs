using System;
using System.Collections.Generic;
using System.Linq;
//using TriangleNet;
//using TriangleNet.Geometry;
//using TriangleNet.Meshing;
//using TriangleNet.Smoothing;
//using TriangleNet.Tools;
using System.IO;


namespace CodendOOP
{
    abstract class Panel
    {
        //=========================
        // variables
        //=========================

        static int nextID = 0;
        public readonly int ID = 0;

        public readonly double LengthInMeshes;
        public readonly double WidthInMeshes;
        public readonly int Orientation;
        public PanelMaterial Material;

        public List<Node> NodeList;
        public List<TriangleElement> TriangleList;
        public List<Node> LeftSelvedgeNodes;
        public List<Node> RightSelvedgeNodes;

        public double Length;
        public double Width;

        private int nAlong;
        private int nAcross;

        //=========================
        // constructor
        //=========================

        public Panel()
        {
            ID = nextID;
            nextID++;
        }

        public Panel(double LengthInMeshes, double WidthInMeshes, int Orientation, PanelMaterial Material)
        {
            ID = nextID;
            nextID++;
            this.LengthInMeshes = LengthInMeshes;
            this.WidthInMeshes = WidthInMeshes;
            this.Orientation = Orientation;
            this.Material = Material;
            CalcUnstretchPanelSize();
        }

        //=========================
        // methods
        //=========================

        private void InitializeFemLists()
        {
            NodeList = new List<Node>();
            TriangleList = new List<TriangleElement>();
            LeftSelvedgeNodes = new List<Node>();
            RightSelvedgeNodes = new List<Node>();
        }

        private void GetAlongAndAcross(MeshSettings meshSettings)
        {
            if (meshSettings.meshingMethod.Equals("ByMesh", StringComparison.InvariantCulture))
            {
                nAlong = (int)Math.Ceiling(LengthInMeshes / meshSettings.MeshPerElemAlong);
                nAcross = (int)Math.Ceiling(WidthInMeshes / meshSettings.MeshPerElemAcross);
            }
            else
            {
                nAlong = meshSettings.ElemAlongPanel;
                nAcross = meshSettings.ElemAcrossPanel;
            }
        }

        public virtual void CalcUnstretchPanelSize()
        {

        }

        public void PrintInfo()
        {
            Console.WriteLine(" * Panel ID: {0}", ID);
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Panel length", LengthInMeshes, "[meshes]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Panel width", WidthInMeshes, "[meshes]");
            Material.PrintInfo();
        }

        //public void StructuredMesh(int nSouth, int nEast, int nNorth, int nWest, double MaxLength, double MaxWidth, QualityOptions quality)
        //{
        //    Length = MaxLength;
        //    Width = MaxWidth;

        //    var poly = new Polygon();

        //    poly.Add(new Contour(new Vertex[4]      // Add the outer box contour with boundary marker ID.
        //    {
        //        new Vertex(0.0, 0.0, 1),
        //        new Vertex(WidthInMeshes, 0.0, 1),
        //        new Vertex(WidthInMeshes, LengthInMeshes, 1),
        //        new Vertex(0.0, LengthInMeshes, 1)
        //    }, 1));

        //    for (int i = 1; i < nSouth; i++)
        //    {
        //        poly.Add(new Vertex(i * WidthInMeshes / nSouth, 0, 1));
        //        //Console.WriteLine("{0,10:F3}{1,10:F3}", i * Geom.WidthInMeshes / nSouth, 0);
        //    }

        //    for (int i = 1; i < nEast; i++)
        //    {
        //        poly.Add(new Vertex(WidthInMeshes, i * LengthInMeshes / nEast, 1));
        //        //Console.WriteLine("{0,10:F3}{1,10:F3}", Geom.WidthInMeshes, i * Geom.LengthInMeshes / nEast);
        //    }

        //    for (int i = 1; i < nNorth; i++)
        //    {
        //        poly.Add(new Vertex(WidthInMeshes * (1 - i / (double)nNorth), LengthInMeshes, 1));
        //        //Console.WriteLine("{0,10:F3}{1,10:F3}", Geom.WidthInMeshes * (1 - i / (double)nNorth), Geom.LengthInMeshes);
        //    }

        //    for (int i = 1; i < nWest; i++)
        //    {
        //        poly.Add(new Vertex(0, LengthInMeshes * (1 - i / (double)nWest), 1));
        //        //Console.WriteLine("{0,10:F3}{1,10:F3}", 0, Geom.LengthInMeshes * (1 - i / (double)nWest));
        //    }

        //    var options = new ConstraintOptions()
        //    {
        //        ConformingDelaunay = true,
        //        SegmentSplitting = 1,
        //    };

        //    var mesh = (Mesh)poly.Triangulate(options, quality);

        //    var smoother = new SimpleSmoother();            // Smooth mesh and re-apply quality options.
        //    smoother.Smooth(mesh, 100);
        //    mesh.Refine(quality);
        //    mesh.Renumber();

        //    //====================================================
        //    //====================================================
        //    // Calculate mesh quality
        //    var statistic = new QualityMeasure();

        //    statistic.Update(mesh);

        //    // Use the minimum triangle area for region refinement
        //    double area = 3 * statistic.AreaMinimum;

        //    foreach (var t in mesh.Triangles)
        //    {
        //        t.Area = area;
        //    }

        //    // Use per triangle area constraint for next refinement
        //    quality.VariableArea = true;

        //    // Refine mesh to meet area constraint.
        //    mesh.Refine(quality);

        //    // Smooth once again.
        //    smoother.Smooth(mesh);
        //    //====================================================
        //    //====================================================

        //    foreach (var v in mesh.Vertices)
        //    {
        //        NodeList.Add(new Node(v.ID,
        //                                   MaxWidth * v.X / WidthInMeshes,
        //                                   MaxLength * v.Y / LengthInMeshes,
        //                                   0));
        //    }

        //    foreach (var t in mesh.Triangles)
        //    {
        //        TriangleList.Add(new DiamondNetTriangle(t.ID,
        //                                              NodeList[t.GetVertexID(0)], NodeList[t.GetVertexID(1)], NodeList[t.GetVertexID(2)],
        //                                              new double[] { t.GetVertex(0).X, t.GetVertex(1).X, t.GetVertex(2).X },
        //                                              new double[] { t.GetVertex(0).Y, t.GetVertex(1).Y, t.GetVertex(2).Y }
        //                                              ) { Label = ID, Material = Material });
        //    }
        //}        

        public void CreateUnitMesh(MeshSettings meshSettings)
        {
            InitializeFemLists();
            GetAlongAndAcross(meshSettings);

            // rectangular grid
            for (int i = 0; i < nAlong + 1; i++) // 
            {
                for (int j = 0; j < nAcross + 1; j++)
                {
                    NodeList.Add(new Node(i * (nAcross + 1) + j, j / (double)nAcross, i / (double)nAlong, 0));

                    // label for last ring of points (for BC)
                    if (i == nAlong)
                    {
                        NodeList.Last().Label = 1;
                    }

                    // label for first ring of points (for BC)
                    if (i == 0)
                    {
                        NodeList.Last().Label = 0;
                    }
                }
                LeftSelvedgeNodes.Add(NodeList[i * (nAcross + 1) + 0]);
                RightSelvedgeNodes.Add(NodeList[i * (nAcross + 1) + nAcross]);
            }

            int gridSize = NodeList.Count;
            int nextrowIdx, thisrowIdx;
            double xMid, yMid;

            // create midpoints in the square cells and define triangles
            for (int i = 0; i < nAlong; i++)
            {
                for (int j = 0; j < nAcross; j++)
                {
                    thisrowIdx = i * (nAcross + 1) + j;
                    nextrowIdx = (i + 1) * (nAcross + 1) + j;

                    // add mid node to a square cell

                    xMid = (NodeList[thisrowIdx].X + NodeList[thisrowIdx + 1].X + NodeList[nextrowIdx].X + NodeList[nextrowIdx + 1].X) / 4;
                    yMid = (NodeList[thisrowIdx].Y + NodeList[thisrowIdx + 1].Y + NodeList[nextrowIdx].Y + NodeList[nextrowIdx + 1].Y) / 4;

                    NodeList.Add(new Node(gridSize + i * nAcross + j, xMid, yMid, 0));

                    // label for pre-last ring of points (for initial shape)
                    if (i == nAlong - 1)
                    {
                        NodeList.Last().Label = 2;
                    }

                    //south trinagle in this row
                    TriangleList.Add(new DiamondNetTriangle(i * 4 * nAcross + 4 * j + 0,
                                    NodeList[gridSize + i * nAcross + j], NodeList[thisrowIdx], NodeList[thisrowIdx + 1],
                                    TriangleuCoord(NodeList, gridSize + i * nAcross + j, thisrowIdx, thisrowIdx + 1),
                                    TrianglevCoord(NodeList, gridSize + i * nAcross + j, thisrowIdx, thisrowIdx + 1))
                    { Material = Material, orientation = Orientation });

                    // east trinagle in this row
                    TriangleList.Add(new DiamondNetTriangle(i * 4 * nAcross + 4 * j + 1,
                                    NodeList[gridSize + i * nAcross + j], NodeList[thisrowIdx + 1], NodeList[nextrowIdx + 1],
                                    TriangleuCoord(NodeList, gridSize + i * nAcross + j, thisrowIdx + 1, nextrowIdx + 1),
                                    TrianglevCoord(NodeList, gridSize + i * nAcross + j, thisrowIdx + 1, nextrowIdx + 1))
                    { Material = Material, orientation = Orientation });

                    // north trinagle in this row
                    TriangleList.Add(new DiamondNetTriangle(i * 4 * nAcross + 4 * j + 2,
                                    NodeList[gridSize + i * nAcross + j], NodeList[nextrowIdx + 1], NodeList[nextrowIdx],
                                    TriangleuCoord(NodeList, gridSize + i * nAcross + j, nextrowIdx + 1, nextrowIdx),
                                    TrianglevCoord(NodeList, gridSize + i * nAcross + j, nextrowIdx + 1, nextrowIdx))
                    { Material = Material, orientation = Orientation });

                    // west trinagle in this row
                    TriangleList.Add(new DiamondNetTriangle(i * 4 * nAcross + 4 * j + 3,
                                    NodeList[gridSize + i * nAcross + j], NodeList[nextrowIdx], NodeList[thisrowIdx],
                                    TriangleuCoord(NodeList, gridSize + i * nAcross + j, nextrowIdx, thisrowIdx),
                                    TrianglevCoord(NodeList, gridSize + i * nAcross + j, nextrowIdx, thisrowIdx))
                    { Material = Material, orientation = Orientation });
                }
            }
        }

        public void LabelNodesForStrap(RoundStrap roundStrap)
        {
            double stepSize = 1.0 / nAcross;
            double roundedPosition = Math.Round( (1 - roundStrap.position) / stepSize) * stepSize;

            foreach (var n in NodeList)
            {
                if (Math.Round(n.Y,2) == Math.Round(roundedPosition,2))
                {
                    n.Label = roundStrap.Label;
                }
            }
        }

        public virtual double[] TriangleuCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { WidthInMeshes * NormalizedXY[idx1].X,
                                  WidthInMeshes * NormalizedXY[idx2].X,
                                  WidthInMeshes * NormalizedXY[idx3].X };
        }

        public virtual double[] TrianglevCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { LengthInMeshes * NormalizedXY[idx1].Y,
                                  LengthInMeshes * NormalizedXY[idx2].Y,
                                  LengthInMeshes * NormalizedXY[idx3].Y };
        }

        /* Mesh operations */

        public void ScaleUnitMesh(double MaxLength, double MaxWidth)
        {
            Width = MaxWidth;
            Length = MaxLength;
            // scale the nodes to deserved distance
            for (int i = 0; i < NodeList.Count; i++)
            {
                NodeList[i].X *= MaxWidth;
                NodeList[i].Y *= MaxLength;
            }

            // scaling changed the geometry, recalculate U V and d of a triangle
            for (int i = 0; i < TriangleList.Count; i++)
            {
                TriangleList[i].UpdateElement();
            }
        }

        public void Translate(double dx, double dy, double dz)
        {
            double X;
            double Y;
            double Z;

            for (int i = 0; i < NodeList.Count; i++)
            {
                X = NodeList[i].X;
                Y = NodeList[i].Y;
                Z = NodeList[i].Z;

                NodeList[i].X = X + dx;
                NodeList[i].Y = Y + dy;
                NodeList[i].Z = Z + dz;
            }
        }

        public void RotateAroundZ(double poleX, double poleZ, double angle)
        {
            double X;
            double Z;

            for (int i = 0; i < NodeList.Count; i++)
            {
                X = NodeList[i].X;
                Z = NodeList[i].Z;

                NodeList[i].X = poleX + (X - poleX) * Math.Cos(angle) - (Z - poleZ) * Math.Sin(angle);


                NodeList[i].Z = poleZ + (X - poleX) * Math.Sin(angle) + (Z - poleZ) * Math.Cos(angle);

            }
        }

        public void MapToCyllinder(double radius, double totalAngle, double startAngle)
        {
            double minX = NodeList[0].X;

            for (int i = 0; i < NodeList.Count; i++)
            {
                if (minX > NodeList[i].X)
                {
                    minX = NodeList[i].X;
                }
            }

            double scale = totalAngle * radius / Width;

            Translate(-minX, 0, 0);

            double theta;

            for (int i = 0; i < NodeList.Count; i++)
            {
                theta = startAngle + NodeList[i].X * scale / radius;
                NodeList[i].X = radius * Math.Cos(theta);
                NodeList[i].Z = radius * Math.Sin(theta);
            }

        }

        public void MapToAxiModel(double[] AxiX, double totalAngle, double startAngle)
        {
            double idxThis = 0;
            int idxPrev, idxNext;
            double R,theta,U, V;

            for (int i = 0; i < NodeList.Count; i++)
            {
                V = NodeList[i].Y * LengthInMeshes;
                U = NodeList[i].X * WidthInMeshes;
               
                idxThis = V * 2;
                idxPrev = Convert.ToInt32(Math.Floor(idxThis));
                idxNext = Convert.ToInt32(Math.Ceiling(idxThis));

                //Console.WriteLine("{0} {1} {2} {3} {4}", NodeList[i].Y, V, idxThis,idxNext,idxPrev);
                NodeList[i].Y = AxiX[2 * idxPrev] + (AxiX[2 * idxNext] - AxiX[2 * idxPrev]) * (idxThis - idxPrev);
                R = AxiX[2 * idxPrev + 1] + (AxiX[2 * idxNext + 1] - AxiX[2 * idxPrev + 1]) * (idxThis - idxPrev);

                theta = startAngle + totalAngle * U / WidthInMeshes;
                NodeList[i].X = R * Math.Cos(theta);
                NodeList[i].Z = R * Math.Sin(theta);
            }
        }

    }
}

