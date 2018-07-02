using System;
using System.Collections.Generic;
using System.IO;

namespace CodendOOP
{
    class RoundStrap
    {
        // instructions at which number of meshes it should be fixed (influences unit panel discretization)
        // instructions how to create bars and whats their L0 based on round strap total length
        //=========================
        // variables
        //=========================

        public int ID;
        public List<Node> NodeList = new List<Node>();
        public List<BarElement> BarList = new List<BarElement>();
        public double Length;
        public double position;
        public BarMaterial Material;
        public int Label {get; set;}
        //=========================
        // constructor
        //=========================

        public RoundStrap(int ID, double position, double Length, BarMaterial Material)
        {
            this.ID = ID;
            this.position = position;
            this.Length = Length;
            this.Material = Material;
        }

        public RoundStrap(int ID, InputReader inputReader)
        {
            this.ID = ID;
            Label = 100 * ID;
            LoadInput(ID, inputReader.Input3d);
        }

        //=========================
        // methods
        //=========================

        public void SetNodeList(List<Node> newList)
        {
            NodeList = new List<Node>(newList);
        }

        public void CreateBarElements()
        {
            InitializeFemLists();

            if (NodeList.Count != 0)
            {
                int BarCount = NodeList.Count - 1;
                double L0 = Length / BarCount;

                for (int i = 0; i < BarCount; i++)
                {
                    BarList.Add(new BarElement(i, NodeList[i], NodeList[i + 1], L0) { Material = Material });
                }
            }
            else
            {
                throw new MissingFieldException("Node list is not set");
            }
        }

        private void InitializeFemLists()
        {
            BarList = new List<BarElement>();
        }

        private void LoadInput(int strapID, string inputPath)
        {
            string targetWord = "RoundStrapCount";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(inputPath);
            string[] parts;
            bool targetFound = false;

            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    targetFound = true;
                    break;
                }
                currentLine++;
            }

            if (targetFound)
            {
                currentLine ++; // skip header
                parts = lines[currentLine + strapID].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                Length = Convert.ToDouble(parts[1]);
                position = Convert.ToDouble(parts[2]);

                if (position < 0 || position >= 1)
                {
                    throw new ArgumentException("Position input for round strap should be between 0 and 1");
                }

                Material = new BarMaterial()
                {
                    Density = Convert.ToDouble(parts[3]),
                    Diameter = Convert.ToDouble(parts[4]),
                    EA = Convert.ToDouble(parts[5]),
                    EI = Convert.ToDouble(parts[6])
                };
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }
    }
}
