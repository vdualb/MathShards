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

using static MathShards.Quadrature.Gauss;

namespace MathShards.FiniteElements.Line.Hermit;

public static class Cubic
{
    public static readonly Func<Real, Real>[] BasisTemplate =
    {
        a => 1 - 3*a*a + 2*a*a*a,
        a => a - 2*a*a + a*a*a,
        a => 3*a*a - 2*a*a*a,
        a => -a*a + a*a*a
    };
    
    public static Real BasisConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisTemplate[i](p01);
    }
    
    public static readonly Func<Real, Real>[] BasisGradTemplate =
    {
        a => -6*a + 6*a*a,
        a => 1 - 4*a + 3*a*a,
        a => 6*a - 6*a*a,
        a => -2*a + 3*a*a
    };
    
    public static Real BasisGradConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisGradTemplate[i](p01);
    }
    
    public static readonly Func<Real, Real>[] BasisGradGradTemplate =
    {
        a => -6 + 12*a,
        a => -4 + 6*a,
        a => 6 - 12*a,
        a => -2 + 6*a
    };
    
    public static Real BasisGradGradConverted(int i, Real p0, Real p1, Real p)
    {
        Real h = p1 - p0;
        Real p01 = (p - p0)/h;
        // см. кирпич с.151
        Real[] coeffs = [1, h, 1, h];

        return coeffs[i] * BasisGradGradTemplate[i](p01);
    }
}
