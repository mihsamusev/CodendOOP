using System;
using System.Collections.Generic;
using static System.Console;

namespace CodendOOP
{
    abstract class NonLinElement
    {
        // Variables
        static int nextID = 0;
        public readonly int ID = 0;
        public int[] globalDOF;
        public int dofPerNode = 3;
        public int nodePerElem;
        public int dofPerElem;
        public List<Node> ElemNodes; 
        // Constructor

        public NonLinElement()
        {
            ID = nextID;
            nextID++;
        }

        public NonLinElement(List<Node> ElemNodes) // constructor without material properties
        {
            ID = nextID;
            nextID++;
            this.ElemNodes = ElemNodes;
            nodePerElem = ElemNodes.Count;
            dofPerElem = 3 * nodePerElem;
            GetGlobalDOF();
        }

        // Methods
        public virtual void PrintInfo()
        {
            WriteLine("Element description goes here");
        }

        public void GetGlobalDOF()
        {
            globalDOF = new int[dofPerElem];

            for (int i = 0; i < nodePerElem; i++)
            {
                for (int j = 0; j < dofPerNode; j++)
                {
                    globalDOF[dofPerNode * i + j] = dofPerNode * ElemNodes[i].ID + j;
                }
            }
        }

        public virtual double[] GetTotalForces()
        {
            return null;
        }

        public virtual double[,] GetTotalStiffness()
        {
            return null;
        }

        public virtual void UpdateElement()
        {

        }
    }
}
