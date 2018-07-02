using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CodendOOP
{
    class BarMaterial
    {
        //==================
        // fields
        //==================

        public double Diameter { get; set; }
        public double EA { get; set; }
        public double EI { get; set; }
        public double Density { get; set; }

        //=================
        // constructor
        //=================

        /*Default*/

        //=================
        // methods
        //=================

        public void PrintInfo()
        {
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Twine thickness", Diameter, "[m]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Density", Density, "[kg / m3]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "EA", EA, "[N]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "EI", EI, "[N * m^2]");
            Console.WriteLine();
        }
    }
}
