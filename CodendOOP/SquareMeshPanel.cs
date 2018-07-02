using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class SquareMeshPanel : Panel
    {
        //=========================
        // variables
        //=========================


        //=========================
        // constructors
        //=========================

        public SquareMeshPanel(double LengthInMeshes, double WidthInMeshes, int Orientation, PanelMaterial Material) : base(LengthInMeshes, WidthInMeshes, Orientation, Material)
        {
            if (!Material.MeshType.Equals("Square",StringComparison.InvariantCultureIgnoreCase))
            {
                throw new ArgumentException("Material mesh type should be square for square mesh panel");
            }
        }

        //=========================
        // methods
        //=========================

        public override void CalcUnstretchPanelSize()
        {
            Length = LengthInMeshes * Material.MeshSide / 2;
            Width = WidthInMeshes * Material.MeshSide / 2;
        }

        public override double[] TriangleuCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { (WidthInMeshes * NormalizedXY[idx1].X - LengthInMeshes * NormalizedXY[idx1].Y) / 2,
                                  (WidthInMeshes * NormalizedXY[idx2].X - LengthInMeshes * NormalizedXY[idx2].Y) / 2,
                                  (WidthInMeshes * NormalizedXY[idx3].X - LengthInMeshes * NormalizedXY[idx3].Y) / 2 };
        }

        public override double[] TrianglevCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { (WidthInMeshes * NormalizedXY[idx1].X + LengthInMeshes * NormalizedXY[idx1].Y) / 2,
                                  (WidthInMeshes * NormalizedXY[idx2].X + LengthInMeshes * NormalizedXY[idx2].Y) / 2,
                                  (WidthInMeshes * NormalizedXY[idx3].X + LengthInMeshes * NormalizedXY[idx3].Y) / 2 };
        }
    }
}
