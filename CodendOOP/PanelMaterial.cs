using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace CodendOOP
{
    class PanelMaterial
    {
        //==================
        // fields
        //==================

        public int ID { get; set; }
        public string MeshType { get; set; }
        public double MeshSide { get; set; }
        public double TwineThickness { get; set; }
        public bool IsDoubleTwine { get; set; }
        public double InitialOpeningAngle { get; set; }
        public double MinimumOpeningAngle { get; set; }
        public double EA { get; set; }
        public double EI { get; set; }
        public double OpenningStifness { get; set; }
        public double KnotContactStifness { get; set; }
        public double Density { get; set; }

        //=================
        // constructor
        //=================

        public PanelMaterial(int ID, InputReader inputReader)
        {
            this.ID = ID;
            LoadInput(ID, inputReader.Input3d);
        }

        //=================
        // methods
        //=================

        private void LoadInput(int materialID, string inpPath)
        {
            string targetWord = "PanelMaterialCount";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            bool targetFound = false;

            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    targetFound = true;
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    if (materialID < 1 || materialID > Convert.ToInt16(parts[1]))
                    {
                    throw new ArgumentException("Wrong panel material ID");
                    }
                    else
                    {
                        break;
                    }
                }
                currentLine++;
            }

            if (targetFound)
            {
                currentLine++; // skip header
                parts = lines[currentLine + materialID].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                MeshType = parts[1];
                Density = Convert.ToDouble(parts[2]);
                MeshSide = Convert.ToDouble(parts[3]);
                IsDoubleTwine = Convert.ToBoolean(Convert.ToInt32(parts[4]));
                TwineThickness = Convert.ToDouble(parts[5]);
                InitialOpeningAngle = Convert.ToDouble(parts[6]);
                EA = Convert.ToDouble(parts[7]);
                EI = Convert.ToDouble(parts[8]);
                OpenningStifness = Convert.ToDouble(parts[9]);
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        public void PrintInfo()
        {
            Console.WriteLine("{0,-25}{1,-10:D}", "Material ID", ID);
            Console.WriteLine("{0,-25}{1,-10:S}", "Mesh type", MeshType);
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Mesh side", MeshSide, "[m]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Twine thickness", TwineThickness, "[m]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Initial opening angle", InitialOpeningAngle, "[deg]");
            if (IsDoubleTwine)
            {
                Console.WriteLine("The mesh is double twine");
            }
            else
            {
                Console.WriteLine("The mesh is single twine");
            }
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Density", Density, "[kg / m3]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "EA", EA,"[N]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "EI", EI,"[N * m^2]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Opening Stiffness", OpenningStifness,"[N * rad]");
            Console.WriteLine();
        }


    }
}
