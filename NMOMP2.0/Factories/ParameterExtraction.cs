using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NMOMP2._0.Interfaces.ParameterExtractor;

namespace NMOMP2._0.Factories
{
    class ParameterExtraction
    {
        static public IParameterExtractor GetExtractor() {
            return new DefaultExtractor();
        }
    }
}
