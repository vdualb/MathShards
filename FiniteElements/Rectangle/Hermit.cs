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

using Real = double;

namespace MathShards.FiniteElements.Rectangle.Hermit;

using Dim1 = FiniteElements.Line.Hermit.Cubic;
using MathShards.TelmaCore;
using static MathShards.Quadrature.Gauss;
using MathShards.CoordSystem.Dim2;
using MathShards.Fem.Common;

public static class Cubic
{
    // 16 функций
    public static Func<PairF64, Real> BasisTemplate(int i)
    {
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;

        return pair
            => Dim1.BasisTemplate[mu](pair.X) * Dim1.BasisTemplate[nu](pair.Y);
    }
    
    public static Real BasisConverted(int i, PairF64 p0, PairF64 p1, PairF64 p)
    {
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;
        
        return Dim1.BasisConverted(mu, p0.X, p1.X, p.X)
            * Dim1.BasisConverted(nu, p0.Y, p1.Y, p.Y);
    }
    
    /// xy==0 => diff x
    /// xy==1 => diff y
    public static Func<PairF64, Real> BasisGradTemplate(int i, int xy)
    {
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;

        if (xy == 0)
        {
            return pair
                => Dim1.BasisGradTemplate[mu](pair.X) * Dim1.BasisTemplate[nu](pair.Y);
        } else if (xy == 1) {
            return pair
                => Dim1.BasisTemplate[mu](pair.X) * Dim1.BasisGradTemplate[nu](pair.Y);
        } else {
            throw new ArgumentException("Invalid xy argument");
        }
    }
    
    /// xy==0 => diff x
    /// xy==1 => diff y
    public static Real BasisGradConverted(int i, int xy, PairF64 p0, PairF64 p1, PairF64 p)
    {
        throw new NotImplementedException("Not finished");
        int mu = 2*(i/4%2) + i%2;
        int nu = 2*(i/8) + i/2%2;

        if (xy == 0)
        {
            return Dim1.BasisGradConverted(mu, p0.X, p1.X, p.X) * Dim1.BasisConverted(mu, p0.Y, p1.Y, p.Y);
        } else if (xy == 1) {
            return Dim1.BasisConverted(mu, p0.X, p1.X, p.X) * Dim1.BasisGradConverted(mu, p0.Y, p1.Y, p.Y);
        } else {
            throw new ArgumentException("Invalid xy argument");
        }
    }
    
    public static Real[,] ComputeLocal<Tc>(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    where Tc : ICoordSystem
    {
        throw new NotImplementedException("I have a great suspicion this doesn't work");
        // side
        var sd = 16;
        var values = new Real[sd, sd];

        for (int i = 0; i < sd; i++)
        {
            for (int j = 0; j < sd; j++)
            {
                var ph = p1 - p0;
                
                var funcMass = (PairF64 point) => {
                    return funcs.Gamma(subDom, point.X, point.Y)
                        * BasisConverted(i, p0, p1, point)
                        * BasisConverted(j, p0, p1, point)
                        * Tc.Jacobian(point.X, point.Y);
                };
                
                values[i, j] = Integrate2DOrder5(p0, p1, funcMass);
                
                var funcStiffness = (PairF64 point) => {
                    // в координатах шаблонного базиса - [0;1]
                    var p01 = new PairF64 (
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    // TODO: не BasisGradTemplate а BasisGradConverted
                    return  funcs.Lambda(subDom, point.X, point.Y)
                        *
                        (
                            BasisGradTemplate(i, 0)(p01)
                            * BasisGradTemplate(j, 0)(p01) / ph.X / ph.X
                        +
                            BasisGradTemplate(i, 1)(p01)
                            * BasisGradTemplate(j, 1)(p01) / ph.Y / ph.Y
                        )
                        * Tc.Jacobian(point.X, point.Y);
                };
                
                values[i, j] += Integrate2DOrder5(p0, p1, funcStiffness);
            }
        }

        return values;
    }
    
    public static Real[] ComputeLocalB<Tc>(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    where Tc : ICoordSystem
    {
        throw new NotImplementedException("I have a great suspicion this doesn't work");
        // side
        var sd = 16;
        var ph = p1 - p0;
        var res = new Real[sd];
        
        for (int i = 0; i < sd; i++)
        {
            var func = (PairF64 point) =>
            {
                return funcs.F(subDom, point.X, point.Y)
                    * BasisConverted(i, p0, p1, point)
                    * Tc.Jacobian(point.X, point.Y);
            };
            res[i] = Integrate2DOrder5(p0, p1, func);
        }
        
        return res;
    }
}
