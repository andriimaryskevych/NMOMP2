using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NMOMP2._0
{
    class Testing
    {
        public Testing()
        {
            // #1: Testing AKT and NT

            //FiniteElementMethod solve = new FiniteElementMethod(4, 4, 4, 4, 4, 4);
            //solve.Start();
            //int[][] nt = solve.NT;
            //double[][] akt = solve.AKT;
            //double[] point;

            //for (int i = 0; i < 20; i++)
            //{
            //    point = akt[nt[0][i]];
            //    Console.WriteLine($"{i}:  {nt[0][i]} --- ({point[0]}, {point[1]}, {point[2]})");
            //}


            // #2: Testing modified SquareGenerator

            //double[][] fi = new SquareGenerator(-1, 0, 1).getItInSingleDArray();
            //Console.WriteLine(fi.Length);
            //for (int i = 0; i < fi.Length; i++)
            //{
            //    Console.WriteLine($"{fi[i][0],2} {fi[i][1],2} {fi[i][2],2}");
            //}


            // #3: Testing Functions

            //double[][] fi = new SquareGenerator(-1, 0, 1).getItInSingleDArray();
            //for (int i = 0; i < 20; i++)
            //{
            //    double[] coord = fi[i];
            //    for (int j = 0; j < 20; j++)
            //    {
            //        Console.Write(FI.getFi(j, coord[0], coord[1], coord[2]));
            //    }
            //    Console.WriteLine();
            //}


            // #4: Gauss nodes

            //double[][] val = Globals.GaussNodes;
            //for (int i = 0; i < 27; i++)
            //{
            //    Console.WriteLine($"{val[i][0],15} {val[i][1],15} {val[i][2],15}");
            //}

            // #5: DFIABG
            //double[,,] dfiabg = Globals.DFIABG;
            //for (int i = 0; i < 27; i++)
            //{
            //    Console.WriteLine(dfiabg[i, 0, 10]);
            //}
        }
    }
}
