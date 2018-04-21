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
    }
}
