using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Web.Script.Serialization;
using Newtonsoft.Json;

namespace NMOMP2._0
{
    class Program
    {
        static void Main(string[] args)
        {
            int x = int.Parse(args[0]);
            int y = int.Parse(args[1]);
            int z = int.Parse(args[2]);

            int elem_x = int.Parse(args[3]);
            int elem_y = int.Parse(args[4]);
            int elem_z = int.Parse(args[5]);

            double puasson = double.Parse(args[6]);
            double jung = double.Parse(args[7]);
            double pressure = double.Parse(args[8]);

            //int x = 1;
            //int y = 1;
            //int z = 1;
            //int size = 3;
            //int elem_x = size;
            //int elem_y = size;
            //int elem_z = size;

            //double puasson = 0.3;
            //double jung = 1;
            //double pressure = -0.3;

            FiniteElementMethod solve = new FiniteElementMethod(x, y, z, elem_x, elem_y, elem_z, puasson, jung, pressure);
            solve.Start();
        }
    }
}
