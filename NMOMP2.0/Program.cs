using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NMOMP2._0
{
    class Program
    {
        static void Main(string[] args)
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

            double[][] fi = new SquareGenerator(-1, 0, 1).getItInSingleDArray();
            Console.WriteLine(fi.Length);
            for (int i = 0; i < fi.Length; i++)
            {
                Console.WriteLine($"{fi[i][0]} {fi[i][1]} {fi[i][2]}");
            }
        }
    }
}
