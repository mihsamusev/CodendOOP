using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class Node : ICloneable
    {
        //variables
        public int ID = 0;
        public double X;
        public double Y;
        public double Z;
        public int Label { get; set; }
        public int CopyOf = -1; // unique node
        //constructor

        public Node(int ID, double X, double Y, double Z)
        {
            this.ID = ID;
            this.X = X;
            this.Y = Y;
            this.Z = Z;
            Label = -1;
        }

        //methods
        public object Clone()
        {
            return MemberwiseClone();
        }

        public void PrintInfo()
        {
            Console.WriteLine("{0,5}{1,10:F3}{2,10:F3}{3,10:F3}{4,10}{5,5}", ID, X, Y, Z, "label",Label);
        }

        public int GetDof(int i)
        {
            if (i >= 0 && i <3)
            {
                return 3 * ID + i;
            }
            return -1;
        }

        public bool IsEqual(Node other)
        {
            if (Math.Abs(X - other.X) < 1e-6 &&
                Math.Abs(Y - other.Y) < 1e-6 &&
                Math.Abs(Z - other.Z) < 1e-6)
            {
                return true;
            }
            return false;
        }
    }
}
