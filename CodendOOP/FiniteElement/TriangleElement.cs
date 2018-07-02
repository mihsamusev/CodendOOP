using System;
using System.Collections.Generic;
using System.IO;
using static System.Math;

namespace CodendOOP
{
    abstract class TriangleElement
    {
        //=========================
        // variables
        //=========================

        public int ID;                      // primary variables    
        public Node n1;
        public Node n2;
        public Node n3;
                                        
        public double[] uCoord;             // mesh and twine coordinates of a triangle
        public double[] vCoord;             // never change unless a panel is remeshed
        public double[] UCoord;
        public double[] VCoord;

        public double[] U = new double[3];  // twine direction vectors
        public double[] V = new double[3];  // can change each iteration because depend on node positions
        public double d;

        public int Label { get; set; }      // secondary variables
        public PanelMaterial Material { get; set; }
        public Towing Towing;
        public bool HasCatch = false;
        public int orientation = 0;
        public int[] globalDOF = new int[9];
        public double[] localF = new double[9];
        public double[,] localK = new double[9, 9];

        //=========================
        // constructors
        //=========================

        public TriangleElement(int ID, Node n1, Node n2, Node n3, double[] uCoord, double[] vCoord)
        {
            this.ID = ID;
            this.n1 = n1;
            this.n2 = n2;
            this.n3 = n3;
            this.uCoord = uCoord;
            this.vCoord = vCoord;

            UCoord = new double[] {uCoord[0] + vCoord[0],
                                   uCoord[1] + vCoord[1],
                                   uCoord[2] + vCoord[2]};

            VCoord = new double[] {vCoord[0] - uCoord[0],
                                   vCoord[1] - uCoord[1],
                                   vCoord[2] - uCoord[2]};
            
            GetUVd();
        }

        //=========================
        // methods
        //=========================

        public void PrintInfo()
        {
            Console.Write("{0,-10:D}", ID);
            Console.Write("{0,-10:D}{1,-10:D}{2,-10:D}",n1.ID, n2.ID, n3.ID);
            Console.Write("{0,-10:F2}{1,-10:F2}", uCoord[0], vCoord[0]);
            Console.Write("{0,-10:F2}{1,-10:F2}", uCoord[1], vCoord[1]);
            Console.Write("{0,-10:F2}{1,-10:F2}", uCoord[2], vCoord[2]);       
            Console.Write("{0}\n",HasCatch);
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

            if (n3.CopyOf != -1)
            {
                n3 = NodeList.Find(node => node.ID == n3.CopyOf);
            }

            globalDOF = GetGlobalDOF();
        }

        public void UpdateElement()
        {
            GetUVd();
        }

        public virtual void GetUVd()
        {

        }

        public int[] GetGlobalDOF()
        {
            return new int[9] {n1.GetDof(0),
                               n1.GetDof(1),
                               n1.GetDof(2),
                               n2.GetDof(0),
                               n2.GetDof(1),
                               n2.GetDof(2),
                               n3.GetDof(0),
                               n3.GetDof(1),
                               n3.GetDof(2)};
        }



        public double CoeffPressureDrag(double d, double un)
        {
            if(un == 0)
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



        public virtual double[] GetNettingWeightForces()
                {
                    return null;
                }

        public virtual double[] GetTwineTensionForces()
        {
            return null;
        }
       
        public virtual double[] GetTwineDragForces()
        {
            return null;
        }

        public virtual double[] GetCatchPressureForces()
        {
            return null;
        }

        public virtual double[] GetMeshOpeningForces()
        {
            return null;
        }

        public virtual double[] GetTwineBendingForces()
        {
            return null;
        }
        
        public virtual double[] GetTotalForces()
        {
            return new double[9];
        }

        public void UpdateElementTotalForce()
        {
            Array.Clear(localF, 0, 9);
            localF = GetTotalForces();
        }



        public virtual double[,] GetTwineTensionStiffness()
        {
            return null;
        }

        public virtual double[,] GetTwineDragStiffness()
        {
            return null;
        }

        public virtual double[,] GetCatchPressureStiffness()
        {
            return null;
        }

        public virtual double[,] GetMeshOpeningStiffness()
        {
            return null;
        }

        public virtual double[,] GetTwineBendingStiffness()
        {
            return null;
        }

        public virtual double[,] GetTotalStiffness()
        {
            return new double[9, 9];
        }

        public void UpdateElementTotalStiffness()
        {
            Array.Clear(localK, 0, 81);
            localK = GetTotalStiffness();
        }
    }
}

