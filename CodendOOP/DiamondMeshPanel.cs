using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class DiamondMeshPanel : Panel
    {
        //=========================
        // variables
        //=========================



        //=========================
        // constructors
        //=========================

        public DiamondMeshPanel(double LengthInMeshes, double WidthInMeshes, int Orientation, PanelMaterial Material) : base(LengthInMeshes, WidthInMeshes, Orientation, Material)
        {
            if (!Material.MeshType.Equals("Diamond", StringComparison.InvariantCultureIgnoreCase))
            {
                throw new ArgumentException("Material mesh type should be diamond for diamond mesh panel");
            }
        }

        //=========================
        // methods
        //=========================

        public override void CalcUnstretchPanelSize()
        {
            Length = LengthInMeshes * Material.MeshSide * Math.Sin(Material.InitialOpeningAngle * Math.PI / 180);
            Width = WidthInMeshes * Material.MeshSide * Math.Sin(Material.InitialOpeningAngle * Math.PI / 180);
        }

        public override double[] TriangleuCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { WidthInMeshes * NormalizedXY[idx1].X,
                                  WidthInMeshes * NormalizedXY[idx2].X,
                                  WidthInMeshes * NormalizedXY[idx3].X };
        }

        public override double[] TrianglevCoord(List<Node> NormalizedXY, int idx1, int idx2, int idx3)
        {
            return new double[] { LengthInMeshes * NormalizedXY[idx1].Y,
                                  LengthInMeshes * NormalizedXY[idx2].Y,
                                  LengthInMeshes * NormalizedXY[idx3].Y };
        }
    }
}
