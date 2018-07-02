using System;
using System.Collections.Generic;
using System.IO;

namespace CodendOOP
{
    class Selvedge
    {
        //=========================
        // variables
        //=========================

        public int ID;       
        public List<Node> NodeList = new List<Node>();
        public List<BarElement> BarList = new List<BarElement>();
        public double Length;
        public BarMaterial Material;
        
        //=========================
        // constructor
        //=========================

        public Selvedge(int ID, double Length, BarMaterial Material)
        {
            this.ID = ID;
            this.Length = Length;
            this.Material = Material;
        }

        public Selvedge(int ID, InputReader inputReader)
        {
            this.ID = ID;
            LoadInput(inputReader.Input3d);
        }

        //=========================
        // methods
        //=========================

        public int PanelID(int i)
        {
            if (i == 0)
            {
                return ID / 10;
            }
            else if(i == 1)
            {
                return ID % 10;
            }
            else
            {
                throw new IndexOutOfRangeException("Input should be 0 or 1");
            }
        }

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

        private void LoadInput(string inputPath)
        {
            string targetWord = "IncludeSelvedges";
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
                currentLine += 2; // skip header
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                Length = Convert.ToDouble(parts[0]);
                Material = new BarMaterial()
                {
                    Density = Convert.ToDouble(parts[1]),
                    Diameter = Convert.ToDouble(parts[2]),
                    EA = Convert.ToDouble(parts[3]),
                    EI = Convert.ToDouble(parts[4])
                };
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        private void InitializeFemLists()
        {
            BarList = new List<BarElement>();
        }
    }
}
