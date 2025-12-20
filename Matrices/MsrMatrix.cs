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
using Real = double;

using System.Collections.Concurrent;

namespace MathShards.Matrices;
using Types;

// like CSR matrix, but diagonal is stored in a separate array
public class MsrMatrix : IHalves
{
    public Real[] Elems = [];
    public Real[] Di = [];
    public int[] Ia = [0];
    public int[] Ja = [];

    public int Size => Di.Length;
    Span<Real> IMatrix.Di => Di;

    public SparkAlgos.Types.IMatrix GetComputeMatrix()
    {
        return new SparkAlgos.Matrices.MsrMatrix(new()
        {
            Elems = Elems,
            Di = Di,
            Ia = Ia,
            Ja = Ja,
        });
    }

    public IEnumerable<Real> FlatNonZero()
    {
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];

            int curr_a = start;

            for (_ = 0; curr_a < stop; curr_a++)
            {
                if (Ja[curr_a] > i)
                {
                    break;
                }
                else
                {
                    if (Elems[curr_a] != 0) yield return Elems[curr_a];
                }
            }
            yield return Di[i];
            for (_ = 0; curr_a < stop; curr_a++)
            {
                if (Elems[curr_a] != 0) yield return Elems[curr_a];
            }
        }
    }

    public unsafe void Mul(ReadOnlySpan<Real> vec, Span<Real> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, Ia.Length - 1);
        fixed (Real* _p_v = vec)
        fixed (Real* _p_res = res)
        {
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = Ia[i];
                    int stop = Ia[i + 1];
                    Real dot = Di[i] * p_v[i];
                    for (int a = start; a < stop; a++)
                    {
                        dot += Elems[a] * p_v[Ja[a]];
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_v = vec)
        fixed(Real* p_res = res)
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = Di[i] * p_v[i];
            for (int a = start; a < stop; a++)
            {
                dot += Elems[a] * p_v[Ja[a]];
            }
            p_res[i] = dot;
        }
#endif
    }

    public unsafe void LMul(ReadOnlySpan<double> vec, Span<double> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, Ia.Length - 1);
        fixed (Real* _p_v = vec)
        fixed (Real* _p_res = res)
        {
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = Ia[i];
                    int stop = Ia[i + 1];
                    Real dot = Di[i] * p_v[i];
                    for (int a = start; a < stop && Ja[a] < i; a++)
                    {
                        dot += Elems[a] * p_v[Ja[a]];
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_v = vec)
        fixed(Real* p_res = res)
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = Di[i] * p_v[i];
            for (int a = start; a < stop && Ja[a] < i; a++)
            {
                dot += Elems[a] * p_v[Ja[a]];
            }
            p_res[i] = dot;
        }
#endif
    }

    public unsafe void UMul(ReadOnlySpan<double> vec, Span<double> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, Ia.Length - 1);
        fixed (Real* _p_v = vec)
        fixed (Real* _p_res = res)
        {
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = Ia[i];
                    int stop = Ia[i + 1];
                    Real dot = Di[i] * p_v[i];
                    for (int a = start; a < stop; a++)
                    {
                        if (Ja[a] > i)
                        {
                            dot += Elems[a] * p_v[Ja[a]];
                        }
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_v = vec)
        fixed(Real* p_res = res)
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = Di[i] * p_v[i];
            for (int a = start; a < stop; a++)
            {
                if (Ja[a] > i)
                {
                    dot += Elems[a] * p_v[Ja[a]];
                }
            }
            p_res[i] = dot;
        }
#endif
    }

    public void InvLMul(Span<double> inOut)
    {
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = 0;
            for (int a = start; a < stop && Ja[a] < i; a++)
            {
                dot += Elems[a] * inOut[Ja[a]];
            }
            inOut[i] -= dot;
            inOut[i] /= Di[i];
        }
    }

    public void InvUMul(Span<double> inOut)
    {
        for (int i = Size - 1; i >= 0; i--)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = 0;
            for (int a = start; a < stop; a++)
            {
                if (Ja[a] > i)
                {
                    dot += Elems[a] * inOut[Ja[a]];
                }
            }
            inOut[i] -= dot;
            inOut[i] /= Di[i];
        }
    }
}
