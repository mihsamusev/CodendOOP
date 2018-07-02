using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class Restraint
    {
        // variables
        public int nodeID;
        public bool xIsFixed;
        public bool yIsFixed;
        public bool zIsFixed;

        //constuctor
        public Restraint(int nodeID, bool xIsFixed, bool yIsFixed, bool zIsFixed)
        {
            this.nodeID = nodeID;
            this.xIsFixed = xIsFixed;
            this.yIsFixed = yIsFixed;
            this.zIsFixed = zIsFixed;
        }

        public void PrintInfo()
        {
            Console.WriteLine("{0,5}{1,10:F3}{2,10:F3}{3,10:F3}", nodeID, xIsFixed, yIsFixed, zIsFixed);
        }
    }
}
