﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NMOMP2._0.DTO
{
    class Parameters
    {
        public Parameters() {}

        public int sizeX { get; set; }
        public int sizeY { get; set; }
        public int sizeZ { get; set; }

        public int xAxisFEMCount { get; set; }
        public int yAxisFEMCount { get; set; }
        public int zAxisFEMCount { get; set; }

        public double puasson { get; set; }
        public double jung { get; set; }
        public double pressure { get; set; }
    }
}
