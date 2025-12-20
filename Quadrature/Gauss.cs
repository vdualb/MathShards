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

namespace MathShards.Quadrature;

using TelmaCore;

public static class Gauss
{
    readonly static Quadrature<PairF64> _cache = Get2DOrder5();
    
    /// p0 - нижний левый угол прямоугольной области
    /// p1 - верхний правый угол
    /// func - на промежутке [-1:1]
    public static double Integrate2DOrder5(
        PairF64 p0, PairF64 p1,
        Func<PairF64, double> func
    ) {
        var quad = _cache;
        var hx = p1.X - p0.X;
        var hy = p1.Y - p0.Y;

        var res = 0.0;
        foreach (var node in quad.Nodes)
        {
            var p = node.Point;
            var w = node.Weight;
            
            var x = new PairF64 (
                hx*(p.X+1.0)/2.0 + p0.X,
                hy*(p.Y+1.0)/2.0 + p0.Y
            );
            
            res += func(x) * w;
        }

        // с Якобианом
        return res * hx*hy / 4.0;
    }
    
    public static double Integrate1DOrder5(
        double p0, double p1,
        Func<double, double> func
    ) {
        var quad = Get1DOrder5();
        var h = p1 - p0;
        
        var res = 0.0;
        foreach (var node in quad.Nodes)
        {
            var p = node.Point;
            var w = node.Weight;

            double x = h * (p + 1.0) / 2.0 + p0;

            res += func(x) * w;
        }

        // с Якобианом
        return res * h / 2.0;
    }
    
    static Quadrature<PairF64> Get2DOrder5()
    {
        var dim1 = Get1DOrder5();
        return Make2D(dim1);
    }
    
    // https://w.wiki/6kqB
    static Quadrature<double> Get1DOrder5()
    {
        double[] points = {
            -Math.Sqrt(0.6),
             0.0,
             Math.Sqrt(0.6)
        };
        double[] weights = [
            5 / 9.0,
            8 / 9.0,
            5 / 9.0
        ];

        var res = new Node<double>[3];

        for (int i = 0; i < 3; i++)
        {
            res[i] = new Node<double>(points[i], weights[i]);
        }

        return new Quadrature<double>(res);
    }
    
    static Quadrature<PairF64> Make2D(Quadrature<double> dim1)
    {
        var len = dim1.Nodes.Length;
        var res = new Node<PairF64>[len * len];

        int i = 0;
        foreach (var node1 in dim1.Nodes)
        {
            foreach (var node2 in dim1.Nodes)
            {
                res[i] = new Node<PairF64> (
                    new(node1.Point, node2.Point),
                    node1.Weight * node2.Weight
                );
                i++;
            }
        }

        return new Quadrature<PairF64>(res);
    }
}
