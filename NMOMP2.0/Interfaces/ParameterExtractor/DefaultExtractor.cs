using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NMOMP2._0.DTO;

namespace NMOMP2._0.Interfaces.ParameterExtractor
{
    class DefaultExtractor : IParameterExtractor
    {
        public DefaultExtractor() { }

        public Parameters extract() {
            Parameters parameters = new Parameters();

            parameters.sizeX = 100;
            parameters.sizeY = 100;
            parameters.sizeZ = 100;

            parameters.xAxisFEMCount = 2;
            parameters.yAxisFEMCount = 2;
            parameters.zAxisFEMCount = 2;

            parameters.puasson = 0.3;
            parameters.jung = 1;
            parameters.pressure = -0.3;

            return parameters;
        }
    }
}
