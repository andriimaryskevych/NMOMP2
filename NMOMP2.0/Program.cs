﻿using System;
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
            //Testing test = new Testing();

            FiniteElementMethod solve = new FiniteElementMethod(4, 4, 4, 4, 4, 4, 0.3);
            solve.Start();

        }
    }
}
