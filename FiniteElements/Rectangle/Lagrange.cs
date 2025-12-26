/*
MathShards
Copyright (C) 2025 Afonin Anton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

namespace MathShards.FiniteElements.Rectangle.Lagrange;

using MathShards.Fem.Common;
using MathShards.TelmaCore;
using static MathShards.Quadrature.Gauss;
using MathShards.CoordSystem.Dim2;

public static class BiLinear
{
    public static readonly Func<PairReal, Real>[] Basis =
    {
        vert => (1 - vert.X) * (1 - vert.Y),
        vert => vert.X       * (1 - vert.Y),
        vert => (1 - vert.X) * vert.Y,
        vert => vert.X       * vert.Y
    };

    public static readonly Func<PairReal, Real>[,] BasisGrad =
    {
        {
            vert => -(1 - vert.Y),
            vert => -(1 - vert.X),
        },
        {
            vert => (1 - vert.Y),
            vert => -vert.X,
        },
        {
            vert => -vert.Y,
            vert => (1 - vert.X),
        },
        {
            vert => vert.Y,
            vert => vert.X
        }
    };
    
    public static Real[,] ComputeLocal<Tc>(ITaskFuncs funcs, PairReal p0, PairReal p1, int subDom)
    where Tc : ICoordSystem
    {
        var values = new Real[4, 4];

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                var ph = p1 - p0;

                Real funcMass(PairReal point)
                {
                    // в координатах шаблонного базиса, [0;1]
                    var p01 = new PairReal(
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return funcs.Gamma(subDom, point.X, point.Y)
                        * Basis[i](p01)
                        * Basis[j](p01)
                        * Tc.Jacobian(point.X, point.Y);
                }

                values[i, j] = Integrate2DOrder5(p0, p1, funcMass);

                Real funcStiffness(PairReal point)
                {
                    // в координатах шаблонного базиса - [0;1]
                    var p01 = new PairReal(
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return funcs.Lambda(subDom, point.X, point.Y)
                        *
                        (
                            BasisGrad[i, 0](p01)
                            * BasisGrad[j, 0](p01) / ph.X / ph.X
                        +
                            BasisGrad[i, 1](p01)
                            * BasisGrad[j, 1](p01) / ph.Y / ph.Y
                        )
                        * Tc.Jacobian(point.X, point.Y);
                }

                values[i, j] += Integrate2DOrder5(p0, p1, funcStiffness);
            }
        }

        return values;
    }
    
    public static Real[] ComputeLocalB<Tc>(ITaskFuncs funcs, PairReal p0, PairReal p1, int subDom)
    where Tc : ICoordSystem
    {
        var ph = p1 - p0;
        var res = new Real[4];
        
        for (int i = 0; i < 4; i++)
        {
            var func = (PairReal point) =>
            {
                // в координатах шаблонного базиса - [0;1]
                var p01 = new PairReal(
                    (point.X - p0.X) / ph.X,
                    (point.Y - p0.Y) / ph.Y
                );
                return funcs.F(subDom, point.X, point.Y)
                    * Basis[i](p01)
                    * Tc.Jacobian(point.X, point.Y);
            };
            res[i] = Integrate2DOrder5(p0, p1, func);
        }
        
        return res;
    }
    
    // для декартовой системы координат
    public static readonly Real[,] LocalG1 = {
        { 2, -2,  1, -1},
        {-2,  2, -1,  1},
        { 1, -1,  2, -2},
        {-1,  1, -2,  2},
    };
    public static readonly Real[,] LocalG2 = {
        { 2,  1, -2, -1},
        { 1,  2, -1, -2},
        {-2, -1,  2,  1},
        {-1, -2,  1,  2},
    };
    public static readonly Real[,] LocalM = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4},
    };
    
    public static Real[] ComputeLocalBTempl(ITaskFuncs funcs, PairReal p0, PairReal p1, int subDom)
    {
        var ph = p1 - p0;
        var res = new Real[4];
        
        /* правая часть */
        Real f1 = funcs.F(subDom, p0.X, p0.Y);
        Real f2 = funcs.F(subDom, p1.X, p0.Y);
        Real f3 = funcs.F(subDom, p0.X, p1.Y);
        Real f4 = funcs.F(subDom, p1.X, p1.Y);

        Real hx = p1.X - p0.X;
        Real hy = p1.Y - p0.Y;

        res[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
        res[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
        res[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
        res[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
    
        return res;
    }

    public static Real [,] ComputeLocalTempl(ITaskFuncs funcs, PairReal p0, PairReal p1, int subDom)
    {
        Real GetGammaAverage()
        {
            Real res = funcs.Gamma(subDom, p0.X, p0.Y)
                     + funcs.Gamma(subDom, p1.X, p0.Y)
                     + funcs.Gamma(subDom, p0.X, p1.Y)
                     + funcs.Gamma(subDom, p1.X, p1.Y);

            return res / 4;
        }

        Real GetLamdaAverage()
        {
            Real res = funcs.Lambda(subDom, p0.X, p0.Y)
                     + funcs.Lambda(subDom, p1.X, p0.Y)
                     + funcs.Lambda(subDom, p0.X, p1.Y)
                     + funcs.Lambda(subDom, p1.X, p1.Y);

            return res / 4;
        }
        
        Real hy = p1.Y - p0.Y;
        Real hx = p1.X - p0.X;
        Real l_avg = GetLamdaAverage();
        Real g_avg = GetGammaAverage();
        
        var values = new Real[4, 4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                values[i, j] = l_avg / 6 * (hy / hx * LocalG1[i, j] + hx / hy * LocalG2[i, j])
                    + g_avg / 36 * hx * hy * LocalM[i, j];
            }
        }

        return values;
    }
}
