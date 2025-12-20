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

namespace MathShards.SlaeSolver;

using static MathShards.SlaeSolver.Shared;
using Matrices.Types;

// МСГ
public class CgmHost : ISlaeSolver
{
    int _maxIter;
    Real _eps;

    int _n = 0; // размерность СЛАУ
    Real[] r = [];
    Real[] di_inv = [];
    Real[] mr = [];
    Real[] az = [];
    Real[] z = [];

    public CgmHost(int maxIter, Real eps)
    {
        _maxIter = maxIter;
        // TODO: уменьшение eps чтобы точность ответа былв сравнима с BicgStab
        _eps = eps / 1e+7;
    }

    public static ISlaeSolver Construct(int maxIter, Real eps)
        => new CgmHost(maxIter, eps);

    // Выделить память для временных массивов
    // n - длина каждого массива
    public void AllocateTemps(int n)
    {
        if (n != _n)
        {
            _n = n;

            r = new Real[_n];
            di_inv = new Real[_n];
            mr = new Real[_n];
            z = new Real[_n];
            az = new Real[_n];
        }
    }

    public (Real discrep, int iter) Solve<T>(T matrix, Span<Real> b, Span<Real> x)
    where T: IMatrix
    {
        AllocateTemps(x.Length);

        var _b = b;

        var r = this.r.AsSpan();
        var di_inv = this.di_inv.AsSpan();
        var mr = this.mr.AsSpan();
        var z = this.z.AsSpan();
        var az = this.az.AsSpan();

        // precond
        matrix.Di.CopyTo(di_inv);
        Rsqrt(di_inv);
        // 1.
        matrix.Mul(x, z);
        _b.CopyTo(r);
        Axpy(-1, z, r);
        // 2.
        r.CopyTo(z);
        Vmul(z, di_inv);

        r.CopyTo(mr);
        Vmul(mr, di_inv);
        var mrr0 = Dot(mr, r);
        
        int iter = 0;
        for (; iter < _maxIter; iter++)
        {
            // 3.    
            z.CopyTo(az);
            matrix.Mul(z, az);

            var azz = Dot(az, z);
            var alpha = mrr0 / azz;
            // 4.
            Axpy(alpha, z, x);
            // 5.
            Axpy(-alpha, az, r);
            // 6.
            r.CopyTo(mr);
            Vmul(mr, di_inv);
            var mrr1 = Dot(mr, r);
            var beta = mrr1/mrr0;
            // 7.
            Scale(beta, z);
            Axpy(1, mr, z);

            mrr0 = mrr1;

            var rr = Dot(r, r);
            var bb = Dot(b, b);
            if (rr / bb < _eps)
            {
                break;
            }
        }

        matrix.Mul(x, z);
        _b.CopyTo(r);
        Axpy(-1, z, r);
        // BLAS.axpy(_x.Length, -1, t, r);
        var rr2 = Dot(r, r);

        return (rr2, iter);
    }
}
