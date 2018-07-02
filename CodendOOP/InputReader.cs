using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Reflection;
using System.IO;

namespace CodendOOP
{
    class InputReader
    {
        //=========================
        // fields
        //=========================

        public string Input3d;
        public string InputAxi;
        public string AxiExe;

        public string OutputAxiShapes;
        public string Output3dShapes;
        public string Output3dResults;

        //=========================
        // constructors
        //=========================

        public InputReader()
        {
            Input3d = AssemblyPath() + "\\input.txt";
            LoadFromInput(Input3d);
        }

        public InputReader(string Input3d)
        {
            this.Input3d = Input3d;
            LoadFromInput(Input3d);
        }

        //=========================
        // methods
        //=========================

        public string AssemblyPath()
        {
            string fullPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().GetName().CodeBase);
            return fullPath.Replace("file:\\", "");
        }

        public PanelInput ReadPanelInput(int panelID)
        {
            string targetWord = "PanelCount";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(Input3d);
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
                    if (panelID < 1 || panelID > Convert.ToInt16(parts[1]))
                    {
                        throw new ArgumentException("Wrong panel ID");
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
                currentLine++;
                parts = lines[currentLine + panelID].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                PanelInput panelInput = new PanelInput()
                {
                    type = parts[1],
                    materialID = Convert.ToInt16(parts[2]),
                    meshesAlong = Convert.ToInt16(parts[3]),
                    meshesAcross = Convert.ToInt16(parts[4]),
                    orientation = Convert.ToInt16(parts[5])
                };

                return panelInput;
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        public double ReadCodendEntranceRadius()
        {
            string targetWord = "EntranceRadius";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(Input3d);
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
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                return Convert.ToDouble(parts[1]);
            }
            else
            {
                throw new IOException("Incomplete input file, entrance radius was not found");
            }            
        }

        private int GetInputCount(string targetWord)
        {
            int count = 0;
            int currentLine = 0;
            string[] lines = File.ReadAllLines(Input3d);
            string[] parts;

            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    count = Convert.ToInt32(parts[1]);
                    break;
                }
                currentLine++;
            }
            return count;
        }

        public int CountPanels()
        {
            string targetWord = "PanelCount";
            int count = GetInputCount(targetWord);
            if (count != 0)
            {
                return count;
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        public int CountPanelMaterials()
        {
            string targetWord = "PanelMaterialCount";
            int count = GetInputCount(targetWord);
            if (count != 0)
            {
                return count;
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        public bool SelvedgesIncluded()
        {
            string targetWord = "IncludeSelvedges";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(Input3d);
            string[] parts;
            bool include = false;

            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    include = Convert.ToBoolean(Convert.ToInt16(parts[1]));
                    break;
                }
                currentLine++;
            }
            return include;
        }

        public bool RoundStrapsIncluded()
        {
            string targetWord = "IncludeRoundStraps";
            int currentLine = 0;
            string[] lines = File.ReadAllLines(Input3d);
            string[] parts;
            bool include = false;

            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    include = Convert.ToBoolean(Convert.ToInt16(parts[1]));
                    break;
                }
                currentLine++;
            }
            return include;
        }

        public int CountRoundStraps()
        {
            string targetWord = "RoundStrapCount";
            int count = GetInputCount(targetWord);
            if (count != 0)
            {
                return count;
            }
            else
            {
                throw new IOException(String.Format("Incomplete input file, case sensitive \'{0}\', \'{1}\' or \'{2}\' tag is not found!",
                                                    targetWord, targetWord.ToUpper(), targetWord.ToLower()));
            }
        }

        public void PrintInfo()
        {
            Console.WriteLine("\nPATHS:");
            Console.WriteLine("Input to the current 3d model:");
            Console.WriteLine("\t" + Input3d);
            Console.WriteLine("Input for pre-calculation with axis-symmetric model:");
            Console.WriteLine("\t" + InputAxi);
            Console.WriteLine("Executable file for axis-symmetric model:");
            Console.WriteLine("\t" + AxiExe);
            Console.WriteLine("Output shape from axis-symmetric model:");
            Console.WriteLine("\t" + OutputAxiShapes);
            Console.WriteLine("Output shape from the current 3d model:");
            Console.WriteLine("\t" + Output3dShapes);
            Console.WriteLine("Output result parameters from the current 3d model:");
            Console.WriteLine("\t" + Output3dResults);
            Console.WriteLine();
        }

        private void LoadFromInput(string input)
        {
            string[] names = { "InputAxi", "AxiExe", "OutputAxiShapes", "Output3dShapes", "Output3dResults" };
            string[] lines = File.ReadAllLines(input);
            string[] parts;
            int currentLine = 0;
            int pathsFound = 0;

            foreach (var line in lines)
            {
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                if (line.Contains(names[0]))
                { InputAxi = parts[1]; pathsFound++; }

                if (line.Contains(names[1]))
                { AxiExe = parts[1]; pathsFound++; }

                if (line.Contains(names[2]))
                { OutputAxiShapes = parts[1]; pathsFound++; }

                if (line.Contains(names[3]))
                { Output3dShapes = parts[1]; pathsFound++; }

                if (line.Contains(names[4]))
                { Output3dResults = parts[1]; pathsFound++; }

                currentLine++;
            }

            if (pathsFound != names.Length)
            {
                throw new Exception("One or more paths could not be found");
            }
        }
    }

    class PanelInput
    {
        public string type;
        public int materialID;
        public int meshesAlong;
        public int meshesAcross;
        public int orientation;
    }
}

