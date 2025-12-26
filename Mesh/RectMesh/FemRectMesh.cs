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

namespace MathShards.Mesh.RectMesh;

using MathShards.Fem.Common;

// TODO: можно ли обойтись без наследования?
public class FemRectMesh : RectMesh
{
    Subdomain[] _subDomains;
    public Subdomain[] SubDomains { get => _subDomains; }

    BoundaryCondition[] _boundaryConditions;
    public BoundaryCondition[] BoundaryConditions { get => _boundaryConditions; }

    public FemRectMesh(
        Real[] xAxis, Real[] yAxis,
        Subdomain[] subDomains,
        BoundaryCondition[] boundaryConditions
    ) : base(xAxis, yAxis) {
        _subDomains = subDomains;
        _boundaryConditions = boundaryConditions;
    }    

    public int? GetSubdomNumAtElCoords (int x1, int y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= IXw[a.X1] && x1 < IXw[a.X2] &&
                y1 >= IYw[a.Y1] && y1 < IYw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }

    public int? GetSubdomNumAtPoint (Real x1, Real y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= Xw[a.X1] && x1 <= Xw[a.X2] &&
                y1 >= Yw[a.Y1] && y1 <= Yw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }
}
