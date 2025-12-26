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

namespace MathShards.SlaeBuilder.Fem;

using MathShards.Fem.Common;
using MathShards.Matrices.Types;
using MathShards.Mesh.RectMesh;

public interface IFemSlaeBuilder
{
    FemRectMesh Mesh { get; }
    GlobalMatrixImplType GlobalMatrixImpl { get; set; }
    static abstract IFemSlaeBuilder Construct(FemRectMesh mesh, ITaskFuncs funcs);
    (IMatrix, Real[]) Build();
}

public enum GlobalMatrixImplType
{
    OpenCL,
    OpenCLV2,
    Host,
    HostParallel,
    HostV2
}

static class Shared {
    
    // https://stackoverflow.com/a/16893641
    public static Real Add(ref Real location1, Real value)
    {
        Real newCurrentValue = location1; // non-volatile read, so may be stale
        while (true)
        {
            Real currentValue = newCurrentValue;
            Real newValue = currentValue + value;
            newCurrentValue = Interlocked.CompareExchange(ref location1, newValue, currentValue);
            if (newCurrentValue.Equals(currentValue))
            {
                return newValue;
            }
        }
    }
    
    // не функция из csharp потому что мне её ещё нужно самому реализовать в OpenCL
    public static int QFind<T> (T[] @where, int start, int end, T what)
    where T: unmanaged, INumber<T>
    {
        int beg = start;
        while (beg < end)
        {
            int mid = (beg + end) / 2;
            if (what > where[mid])
            {
                beg = mid + 1;
            }
            else
            {
                end = mid;
            }
        }

        if (where[beg] != what)
        {
            throw new Exception("Quick search failed");
        }

        return beg;
    }
    
    public static int LFind<T> (T[] where, T what, int start)
    where T: unmanaged, INumber<T>
    {
        while (@where[start] != what) start++;
        return start;
    }

    public static int? TryLFind<T> (T[] where, T what, int start, int end)
    where T: unmanaged, INumber<T>
    {
        while (start < end && @where[start] != what) start++;

        if (start < end)
        {
            return start;
        }
        else
        {
            return null;
        }

    }
}
