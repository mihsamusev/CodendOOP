using System;
using System.IO;

namespace CodendOOP
{
    class MeshSettings
    {
        //============================
        // fields
        //============================

        public int ElemAcrossPanel { get; set; }
        public int ElemAlongPanel { get; set; }
        public double MeshPerElemAcross{ get; set; }
        public double MeshPerElemAlong { get; set; }
        public string meshingMethod { get; set; }

        //============================
        // constructor
        //============================

        public MeshSettings()
        {
            LoadDefault();
        }

        public MeshSettings(InputReader filePath)
        {
            LoadDefault();
            LoadInput(filePath.Input3d);
        }

        //============================
        // methods
        //============================

        private void LoadDefault()
        {
            meshingMethod = "ByElement";
            ElemAcrossPanel = 10;
            ElemAlongPanel = 10;
        }

        private void LoadInput(string inpPath)
        {
            string targetWord = "MeshingMethod";
            string method = "";
            string[] names = { "MeshPerElemAcross", "MeshPerElemAlong", "ElemAcrossPanel", "ElemAlongPanel" };
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            int currentLine = 0;
            bool methodFound = false;


            foreach (var line in lines)
            {
                if (line.Contains(targetWord) ||
                    line.Contains(targetWord.ToUpper()) ||
                    line.Contains(targetWord.ToLower()))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    method = parts[1];
                    methodFound = true;
                    break;
                }
                currentLine++;
            }

            if (!methodFound)
            {
                return;
            }
            else if (methodFound && method.Equals("ByElement", StringComparison.InvariantCultureIgnoreCase))
            {
                meshingMethod = method;
                currentLine = 0;
                foreach (var line in lines)
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                    if (line.Contains(names[2]))
                        ElemAcrossPanel = Convert.ToInt32(parts[1]);

                    if (line.Contains(names[3]))
                        ElemAlongPanel = Convert.ToInt32(parts[1]);

                    currentLine++;
                }
            }
            else if (methodFound && method.Equals("ByMesh", StringComparison.InvariantCultureIgnoreCase))
            {
                meshingMethod = method;
                currentLine = 0;
                foreach (var line in lines)
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                    if (line.Contains(names[0]))
                        MeshPerElemAcross = Convert.ToDouble(parts[1]);

                    if (line.Contains(names[1]))
                        MeshPerElemAlong = Convert.ToDouble(parts[1]);

                    currentLine++;
                }
            }
        }


    }
}
