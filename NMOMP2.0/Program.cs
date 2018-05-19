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
            //Testing test = new Testing();

            FiniteElementMethod solve = new FiniteElementMethod(100, 100, 100, 5, 5, 5, 0.3);
            solve.Start();

            //var v = new { Amount = 108, Message = "Hello" };
            //var json = JavaScriptSerializer.Serialize();

            //double[][] arr = new double[][] { new double[] { 1,2,3}, new double[] { 1, 2, 3 } };

            //string jsonData = JsonConvert.SerializeObject((from a in arr select new { x = a[0], y = a[1], z = a[2], }));

            //Console.WriteLine(jsonData);

        }
    }
}
