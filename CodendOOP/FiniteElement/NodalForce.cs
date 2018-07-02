using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class NodalForce
    {
        // variables
        public int nodeID;
        public double Fx;
        public double Fy;
        public double Fz;

        //constuctor
        public NodalForce(int nodeID, double Fx, double Fy, double Fz)
        {
            this.nodeID = nodeID;
            this.Fx = Fx;
            this.Fy = Fy;
            this.Fz = Fz;
        }

    }
}
