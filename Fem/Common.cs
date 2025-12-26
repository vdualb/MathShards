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

namespace MathShards.Fem.Common;

public struct Subdomain
{
    // номера подобласти начинаются с 0
    public int Num;
    public int X1;
    public int X2;
    public int Y1;
    public int Y2;
}

public struct BoundaryCondition
{
    // в файлах нумерация с 1
    // в программе - с 0
    public int Num { get; set; }
    // тип краевого условия (первый, второй ...) нумеруется с 1
    public int Type { get; set; }
    public int X1 { get; set; }
    public int X2 { get; set; }
    public int Y1 { get; set; }
    public int Y2 { get; set; }
}

public interface ITaskFuncs
{
    string Description { get; }

    Real Answer(int subdom, Real x, Real y);
    Real Lambda(int subdom, Real x, Real y);
    Real Gamma(int subdom, Real x, Real y);
    Real F(int subdom, Real x, Real y);
    
    // I
    Real Ug(int bcNum, Real x, Real y);
    // II
    Real Theta(int bcNum, Real x, Real y);
    // III
    Real Beta(int bcNum);
    Real uBeta(int bcNum, Real x, Real y);
}
