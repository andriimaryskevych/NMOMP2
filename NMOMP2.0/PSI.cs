using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NMOMP2._0
{
    class PSI
    {
        private static double HALF = 0.5;
        private static double ONE_FOURTH = 0.25;

        private static Dictionary<int, double[]> adapter = new Dictionary<int, double[]>
        {
            {0, new double[] {-1, -1} },
            {1, new double[] { 1, -1} },
            {2, new double[] { 1,  1} },
            {3, new double[] {-1,  1} },

            {4, new double[] { 0, -1} },
            {5, new double[] { 1,  0} },
            {6, new double[] { 0,  1} },
            {7, new double[] {-1,  0} }
        };

        private static double diEtaFirst(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = ONE_FOURTH *
                     (1 + tau * coord[1]) *
                     coord[0] *
                     (2 * eta * coord[0] + tau * coord[1]);
            return result;
        }
        private static double diEtaSecond(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = HALF *
                     (1 + tau * coord[1]) *
                     (-2 * eta);
            return result;
        }
        private static double diEtaThird(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = HALF *
                     (1 - Math.Pow(tau, 2)) *
                     coord[0];
            return result;
        }

        private static double diTauFirst(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = ONE_FOURTH *
                     (1 + eta * coord[0]) *
                     coord[1] *
                     (2 * tau * coord[1] + eta * coord[0]);

            return result;
        }
        private static double diTauSecond(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = HALF *
                     (1 - Math.Pow(eta, 2)) *
                     coord[1];
            return result;
        }
        private static double diTauThird(int i, double eta, double tau)
        {
            double result;
            double[] coord = adapter[i];
            result = HALF *
                     (1 + eta * coord[0]) *
                     (-2 * tau);
            return result;
        }

        // functinons delegate calculation to more specified functions depending on their index 
        private static double diEta(int i, double eta, double tau)
        {
            return i < 4 ? diEtaFirst(i, eta, tau) : i % 2 == 0 ? diEtaSecond(i, eta, tau) : diEtaThird(i, eta, tau);
        }
        private static double diTau(int i, double eta, double tau)
        {
            return i < 4 ? diTauFirst(i, eta, tau) : i % 2 == 0 ? diTauSecond(i, eta, tau) : diTauThird(i, eta, tau);
        }

        public static double getDiPsi(int variable, int i, double eta, double tau)
        {
            double result = 0;

            switch (variable)
            {
                case 0: result = diEta(i, eta, tau); break;
                case 1: result = diTau(i, eta, tau); break;
            }

            return result;
        }
    }
}
