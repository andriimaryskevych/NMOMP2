using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace NMOMP2._0
{
    class FiniteElementMethod
    {
        // this one holds all the points
        // when i devide a kube in 2x2x2 parts it's size beacames 5x5x5 * 3
        // to use only primitives, i decidet to store 5x5x5 vectors of length 3
        // this vectors hold exact position of each vertex in 3d space
        public double[,,][] matrix;

        // in case above, lists Count is equal 8
        //public List<FiniteElement> finiteElement;

        public int x;
        public int y;
        public int z;

        public int m;
        public int n;
        public int k;

        public double E = 1;
        public double v;
        public double lam;
        public double mu;


        public int nqp;
        // I changed order here:
        // As first argument I set global number and as a second one i set vector of lenght 3 that contains exact point coord
        public double[][] AKT;

        private int[,,] local_to_global;
        private Dictionary<int, int[]> adapter = Globals.magicDictionary;

        private double[][] GaussNodes = Globals.GaussNodes;
        private double[,,] DFIABG = Globals.DFIABG;

        int nel;
        // Same as AKT, FE number goes first than it's local value that evaluates to global coord of selected FE
        public int[][] NT;

        private double SCALE_X;
        private double SCALE_Y;
        private double SCALE_Z;

        public FiniteElementMethod(int _x, int _y, int _z, int _m, int _n, int _k, double _v)
        {
            x = _x;
            y = _y;
            z = _z;

            m = _m;
            n = _n;
            k = _k;

            v = _v;
            lam = E / ((1 + v) * (1 - 2 * v));
            mu = E / (2 * (1 + v));

            SCALE_X = (double)x / m;
            SCALE_Y = (double)y / n;
            SCALE_Z = (double)z / k;

            matrix = new double[m * 2 + 1, n * 2 + 1, k * 2 + 1][];
            nqp = 0;
            // This is private helper matrix
            // it's purpose is to store global coord related to matrix sizes coord (I can't explain better, just look at createAKT))) )
            local_to_global = new int[m * 2 + 1, n * 2 + 1, k * 2 + 1];

            nel = m * n * k;
            NT = new int[nel][];
        }

        public void Start()
        {
            fillMatrixWithMainVertexes();
            fillMatrixWithIntermidiateVertexes();
            createAKT();
            createNT();
            getMGE(0);
        }

        private void fillMatrixWithMainVertexes()
        {
            for (int deltaZ = 0; deltaZ < k + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m + 1; deltaX++)
                    {
                        matrix[deltaX * 2, deltaY * 2, deltaZ * 2] = new double[] { SCALE_X * deltaX, SCALE_Y * deltaY, SCALE_Z * deltaZ };
                        nqp++;
                    }
                }
            }
        }
        private void fillMatrixWithIntermidiateVertexes()
        {
            double[] current;
            for (int deltaZ = 0; deltaZ < k + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m + 1; deltaX++)
                    {
                        current = matrix[deltaX * 2, deltaY * 2, deltaZ * 2];
                        if (deltaX != m)
                        {
                            matrix[deltaX * 2 + 1, deltaY * 2, deltaZ * 2] = new double[] { current[0] + SCALE_X / 2, current[1], current[2] };
                            nqp++;
                        }

                        if (deltaY != n)
                        {
                            matrix[deltaX * 2, deltaY * 2 + 1, deltaZ * 2] = new double[] { current[0], current[1] + SCALE_Y / 2, current[2] } ;
                            nqp++;
                        }

                        if (deltaZ != k)
                        {
                            matrix[deltaX * 2, deltaY * 2, deltaZ * 2 + 1] = new double[] { current[0], current[1], current[2] + SCALE_Z / 2 };
                            nqp++;
                        }
                    }
                }
            }
        }
        private void createAKT()
        {
            AKT = new double[nqp][];
            int counter = 0;
            for (int deltaZ = 0; deltaZ < k * 2 + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n * 2 + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m * 2 + 1; deltaX++)
                    {
                        if (matrix[deltaX, deltaY, deltaZ] != null)
                        {
                            AKT[counter] = matrix[deltaX, deltaY, deltaZ];
                            local_to_global[deltaX, deltaY, deltaZ] = counter;
                            counter++;
                        }
                    }
                }
            }
        }

        private void createNT()
        {            
            for (int figure = 0; figure < nel; figure++)
            {
                int current_element = figure;
                int Z_COOORDINATE = (int)(figure / (m * n));
                current_element %= (m * n);
                int Y_COORDINATE = (int)(current_element / m);
                current_element %= (m);
                int X_COORDINATE = current_element;

                createIndependentElements(X_COORDINATE * 2, Y_COORDINATE * 2, Z_COOORDINATE * 2, figure);
            }
        }
        private void createIndependentElements(int x, int y, int z, int belongElementNumber)
        {
            int[] globalCoordinates = new int[20];
            int[] delta;
            for (int i = 0; i < 20; i++)
            {
                delta = adapter[i];
                globalCoordinates[i] = local_to_global[x + delta[0], y + delta[1], z + delta[2]];
            }
            NT[belongElementNumber] = globalCoordinates;
        }

        private void getMGE(int number)
        {
            // defined once and used for each finite element
            double[,,] dxyzabg = new double[3, 3, 27];
            double[] dj = new double[27];
            double[,,] dfixyz = new double[27, 20, 3];

            int[] coordinates = NT[number];

            // calc dxyzabg
            double globalCoordinate = 0;
            double diFi = 0;
            double sum = 0;
            
            // i stands for global coordinate
            // j for local
            // k for gaussNode
            // l for i-th function
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 27; k++)
                    {
                        sum = 0;
                        for (int l = 0; l < 20; l++)
                        {
                            globalCoordinate = AKT[coordinates[l]][i];
                            diFi = DFIABG[k, j, l];
                            sum += globalCoordinate * diFi;
                        }
                        dxyzabg[i, j, k] = sum;
                    }
                }
            }

            // calc dj
            double[,] jak;
            for (int i = 0; i < 27; i++)
            {
                jak = new double[3, 3] {
                    { dxyzabg[0,0,i], dxyzabg[1,0,i], dxyzabg[2,0,i] },
                    { dxyzabg[0,1,i], dxyzabg[1,1,i], dxyzabg[2,1,i] },
                    { dxyzabg[0,2,i], dxyzabg[1,2,i], dxyzabg[2,2,i] }
                };
                dj[i] = (
                                jak[0, 0] * jak[1, 1] * jak[2, 2] +
                                jak[0, 1] * jak[1, 2] * jak[2, 0] +
                                jak[0, 2] * jak[1, 0] * jak[2, 1]
                            ) -
                            (
                                jak[0, 2] * jak[1, 1] * jak[2, 0] +
                                jak[0, 1] * jak[1, 0] * jak[2, 2] +
                                jak[0, 0] * jak[1, 2] * jak[2, 1]
                            );
            }

            // calc dfixyz
            // col is free column
            double[] col = new double[3];
            // i stands for gausNode
            // phi stands for i-th function
            for (int i = 0; i < 27; i++)
            {
                for (int phi = 0; phi < 20; phi++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        col[k] = DFIABG[i, k, phi];
                    }

                    double[] gaussianSolve = Gaussian.Solve(new double[3, 3] {
                            { dxyzabg[0,0,i], dxyzabg[1,0,i], dxyzabg[2,0,i] },
                            { dxyzabg[0,1,i], dxyzabg[1,1,i], dxyzabg[2,1,i] },
                            { dxyzabg[0,2,i], dxyzabg[1,2,i], dxyzabg[2,2,i] }
                        }, col);

                    for (int k = 0; k < 3; k++)
                    {
                        dfixyz[i, phi, k] = gaussianSolve[k];
                    }
                }
            }

            // calc mge
            double[,] mge = new double[60, 60];
            for (int i = 0; i < 60; i++)
            {
                for (int j = 0; j < 60; j++)
                {
                    if (i > j)
                    {
                        mge[i, j] = mge[j, i];
                    }
                    else {

                    }
                }
            }
        }

        public double[,] one_one(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        res[i, j] = 11;
                    }
                }
            }
            return res;
        }
        public double[,] two_two(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        res[i, j] = 22;
                    }
                }
            }
            return res;
        }
        public double[,] three_three(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        res[i, j] = 33;
                    }
                }
            }
            return res;
        }

        public double[,] one_two(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    res[i, j] = 12;                    
                }
            }
            return res;
        }
        public double[,] one_three(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    res[i, j] = 13;
                }
            }
            return res;
        }
        public double[,] two_three(double[,,] dfixyz)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    res[i, j] = 23;
                }
            }
            return res;
        }

        public double[,] rotate(double[,] toRotate)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    res[i, j] = toRotate[j, i];
                }
            }
            return res;
        }
    }
}
