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

using System.Numerics;

namespace MathShards.SlaeBuilder.Spline;

using MathShards.Mesh.RectMesh;
using Matrices.Types;

public interface ISplineSlaeBuilder
{
    // RectMesh Mesh { get; }

    /// <param name="inMesh">Сетка из решённой задачи</param>
    /// <param name="inValues">Результат задачи</param>
    /// <param name="mesh">Сетка для построения сплайна</param>
    static abstract ISplineSlaeBuilder Construct(
        RectMesh inMesh, Real[] inValues,
        RectMesh mesh
    );
    (IMatrix matrix, Real[] right) Build();
}

public struct SplineParams
{
    public Real Alpha { get; set; }
    public Real Beta { get; set; }
    public Real W { get; set; }
}

public interface ISplineSlaeBuilder1D
{
    // RectMesh Mesh { get; }

    /// <param name="inMesh">Сетка из решённой задачи</param>
    /// <param name="inValues">Результат задачи</param>
    /// <param name="mesh">Сетка для построения сплайна</param>
    static abstract ISplineSlaeBuilder1D Construct(
        LineMesh inMesh, Real[] inValues,
        LineMesh mesh, SplineParams splineParams
    );
    (IMatrix matrix, Real[] right) Build();
}
