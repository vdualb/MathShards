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

#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

namespace MathShards.Quadrature;

using TelmaCore;

public static class Gauss
{
    readonly static Quadrature<PairReal> _cache = Get2DOrder5();
    
    /// p0 - нижний левый угол прямоугольной области
    /// p1 - верхний правый угол
    /// func - на промежутке [-1:1]
    public static Real Integrate2DOrder5(
        PairReal p0, PairReal p1,
        Func<PairReal, Real> func
    ) {
        var quad = _cache;
        var hx = p1.X - p0.X;
        var hy = p1.Y - p0.Y;

        var res = 0.0;
        foreach (var node in quad.Nodes)
        {
            var p = node.Point;
            var w = node.Weight;
            
            var x = new PairReal (
                (Real)(hx *(p.X+1.0)/2.0 + p0.X),
                (Real)(hy *(p.Y+1.0)/2.0 + p0.Y)
            );
            
            res += func(x) * w;
        }

        // с Якобианом
        return (Real)(res * hx*hy / 4.0);
    }
    
    public static Real Integrate1DOrder5(
        Real p0, Real p1,
        Func<Real, Real> func
    ) {
        var quad = Get1DOrder5();
        var h = p1 - p0;
        
        var res = 0.0;
        foreach (var node in quad.Nodes)
        {
            var p = node.Point;
            var w = node.Weight;

            Real x = (Real)(h * (p + 1.0) / 2.0 + p0);

            res += func(x) * w;
        }

        // с Якобианом
        return (Real)(res * h / 2.0);
    }
    
    static Quadrature<PairReal> Get2DOrder5()
    {
        var dim1 = Get1DOrder5();
        return Make2D(dim1);
    }
    
    // https://w.wiki/6kqB
    static Quadrature<Real> Get1DOrder5()
    {
        Real[] points = [
            -(Real)Math.Sqrt(0.6),
            0.0f,
            (Real)Math.Sqrt(0.6)
        ];
        Real[] weights = [
            (Real)(5 / 9.0),
            (Real)(8 / 9.0),
            (Real)(5 / 9.0)
        ];

        var res = new Node<Real>[3];

        for (int i = 0; i < 3; i++)
        {
            res[i] = new Node<Real>(points[i], weights[i]);
        }

        return new Quadrature<Real>(res);
    }
    
    static Quadrature<PairReal> Make2D(Quadrature<Real> dim1)
    {
        var len = dim1.Nodes.Length;
        var res = new Node<PairReal>[len * len];

        int i = 0;
        foreach (var node1 in dim1.Nodes)
        {
            foreach (var node2 in dim1.Nodes)
            {
                res[i] = new Node<PairReal> (
                    new(node1.Point, node2.Point),
                    node1.Weight * node2.Weight
                );
                i++;
            }
        }

        return new Quadrature<PairReal>(res);
    }
}
