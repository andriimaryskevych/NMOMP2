using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NMOMP2._0
{
    class Globals
    {
        public static Dictionary<int, int[]> magicDictionary = new Dictionary<int, int[]> {
            {0, new int[]{0,0,0} },
            {1, new int[]{2,0,0} },
            {2, new int[]{2,2,0} },
            {3, new int[]{0,2,0} },
            
            {4, new int[]{0,0,2} },
            {5, new int[]{2,0,2} },
            {6, new int[]{2,2,2} },
            {7, new int[]{0,2,2} },

            {8,  new int[]{1,0,0} },
            {9, new int[]{2,1,0} },
            {10, new int[]{1,2,0} },
            {11, new int[]{0,1,0} },
            
            {12, new int[]{0,0,1} },
            {13, new int[]{2,0,1} },
            {14, new int[]{2,2,1} },
            {15, new int[]{0,2,1} },
            
            {16, new int[]{1,0,2} },
            {17, new int[]{2,1,2} },
            {18, new int[]{1,2,2} },
            {19, new int[]{0,1,2} },
        };

        private static double[][] gaussNodes()
        {
            double[][] res = new double[27][];
            double[] val = new double[] { -1 * Math.Sqrt(0.6), 0, Math.Sqrt(0.6) };

            int counter = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        res[counter] = new double[] { val[i], val[j], val[k] };
                        ++counter;
                    }
                }
            }

            return res;
        }
        public static double[][] GaussNodes = gaussNodes();

        private static double[,,] countDFIABG()
        {
            double[,,] result = new double[27, 3, 20];

            // i stands for number of gauss node
            // local stands for local variable, like 0 is alpha...
            // k stands for i-th function
            for (int i = 0; i < 27; i++)
            {
                for (int local = 0; local < 3; local++)
                {
                    for (int k = 0; k < 20; k++)
                    {
                        result[i, local, k] = FI.getDiFi(local, k, GaussNodes[i][0], GaussNodes[i][1], GaussNodes[i][2]);
                    }
                }
            }

            return result;
        }
        public static double[,,] DFIABG = countDFIABG();
    }
}
