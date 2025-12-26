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

#define HOST_PARALLEL
#define USE_BLAS

#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Collections.Concurrent;

using Quasar.Native;

namespace MathShards.SlaeSolver;

public interface ISlaeSolver
{
    static abstract ISlaeSolver Construct(int maxIter, Real eps);
    // x используется как начальное приближение, туда же попадёт ответ
    (Real discrep, int iter) Solve<T> (T matrix, Span<Real> b, Span<Real> x)
    where T: Matrices.Types.IMatrix;
    void AllocateTemps(int n);
}

static class Shared
{
    // y *= x
    public static unsafe void Vmul(Span<Real> y, ReadOnlySpan<Real> x)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, y.Length);
        fixed(Real* _p_y = y)
        fixed(Real* _p_x = x)
        {
            var p_y = _p_y;
            var p_x = _p_x;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    p_y[i] *= p_x[i];
                }
            });
    
        }
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] *= x[i];
        }
#endif
    }
    
    // y = y*(-1/2)
    public static unsafe void Rsqrt(Span<Real> y)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, y.Length);
        fixed(Real* _p_y = y)
        {
            var p_y = _p_y;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    p_y[i] = (Real)(1 / Math.Sqrt(p_y[i]));
                }
            });
    
        }
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] = (Real)(1 / Math.Sqrt(y[i]));
        }
#endif
    }
    
    // y += alpha*x
    public static void Axpy(Real alpha, ReadOnlySpan<Real> x, Span<Real> y)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
#if USE_BLAS
        BLAS.axpy(x.Length, alpha, x, y);
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] += (Real)(alpha * x[i]);
        }
#endif
    }
    // x·y
    public static Real Dot(ReadOnlySpan<Real> x, ReadOnlySpan<Real> y)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
#if USE_BLAS
        return (Real)BLAS.dot(x.Length, x, y);
#else
        Real sum = 0;
        for (int i = 0; i < y.Length; i++)
        {
            sum += x[i] * y[i];
        }
        return sum;
#endif
    }
    // y_i = alpha * y[i]
    public static void Scale(Real alpha, Span<Real> y)
    {
#if USE_BLAS
        BLAS.scal(y.Length, alpha, y);
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] *= alpha;
        }
#endif
    }
}
