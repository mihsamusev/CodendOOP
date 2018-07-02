using System;
using System.IO;

namespace CodendOOP
{
    class Towing
    {
        //=========================
        // variables
        //=========================

        public double Cx;
        public double Cy;
        public double Cz;
        public double[] Vector;
        public double Speed;
        public bool IncludeNettingDrag { get; set; }

        //=========================
        // constructor
        //=========================

        public Towing(double Cx, double Cy, double Cz)
        {
            this.Cx = Cx;
            this.Cy = Cy;
            this.Cz = Cz;
            IncludeNettingDrag = true;
            Update();
        }

        public Towing(InputReader filePath)
        {
            LoadInput(filePath.Input3d);
            Update();
        }

        //=========================
        // methods
        //=========================

        private void LoadInput(string inpPath)
        {
            string[] names = { "CX", "CY", "CZ", "IncludeNettingDrag" };
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            int currentLine = 0;

            foreach (var line in lines)
            {
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                if (line.Contains(names[0]))
                    Cx = Convert.ToDouble(parts[1]);

                if (line.Contains(names[1]))
                    Cy = Convert.ToDouble(parts[1]);

                if (line.Contains(names[2]))
                    Cz = Convert.ToDouble(parts[1]);

                if (line.Contains(names[3]))
                    IncludeNettingDrag = Convert.ToBoolean(Convert.ToInt32(parts[1]));

                currentLine++;
            }
        }

        public void Update()
        {
            Vector = new double[] { Cx, Cy, Cz };
            Speed = Math.Sqrt(Math.Pow(Cx, 2) + Math.Pow(Cy, 2) + Math.Pow(Cz, 2));
        }

    }
}
