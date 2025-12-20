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
public class CgmEisenstatSimpleHost : ISlaeSolver
{
    int _maxIter;
    Real _eps;

    int _n = 0; // размерность СЛАУ
    Real[] r_hat = [];
    Real[] r_stroke = [];
    Real[] p = [];
    Real[] t = [];
    Real[] Ap = [];
    Real[] di_inv = [];
    Real[] mr = [];
    Real[] az = [];
    Real[] z = [];

    public CgmEisenstatSimpleHost(int maxIter, Real eps)
    {
        _maxIter = maxIter;
        // TODO: уменьшение eps чтобы точность ответа былв сравнима с BicgStab
        _eps = eps / 1e+7;
    }

    public static ISlaeSolver Construct(int maxIter, Real eps)
        => new CgmEisenstatSimpleHost(maxIter, eps);

    // Выделить память для временных массивов
    // n - длина каждого массива
    public void AllocateTemps(int n)
    {
        if (n != _n)
        {
            _n = n;

            r_hat =      new Real[_n];
            r_stroke =   new Real[_n];
            p =          new Real[_n];
            t =          new Real[_n];
            Ap =         new Real[_n];
            di_inv =     new Real[_n];
            mr =         new Real[_n];
            az =         new Real[_n];
            z =          new Real[_n];
        }
    }

    
    // https://epubs.siam.org/doi/10.1137/0902001
    public (Real discrep, int iter) Solve<T>(T matrix, Span<Real> b, Span<Real> x)
    where T : IMatrix
    {
        if (matrix is IHalves m)
        {
            return SolveImpl(m, b, x);
        } else {
            throw new ArgumentException();
        }
    }
    
    public (Real discrep, int iter) SolveImpl(IHalves matrix, Span<Real> b, Span<Real> x)
    {
        AllocateTemps(x.Length);

        var r_hat = this.r_hat.AsSpan();
        var r_stroke = this.r_stroke.AsSpan();
        var p = this.p.AsSpan();
        var di_inv = this.di_inv.AsSpan();
        var mr = this.mr.AsSpan();
        var Ap = this.Ap.AsSpan();
        var z = this.z.AsSpan();
        var t = this.t.AsSpan();
        var az = this.az.AsSpan();

        // 2a
        // t_hat but as a part of A_hat*x_0
        matrix.Mul(x, r_stroke);
        b.CopyTo(r_hat);
        Axpy(-1, r_stroke, r_hat);

        // 2b
        r_hat.CopyTo(r_stroke);
        matrix.InvLMul(r_stroke);
        Vmul(r_stroke, matrix.Di);
        matrix.InvUMul(r_stroke);
        
        r_stroke.CopyTo(p);
        
        // precompute rr0
        var rr0 = Dot(r_hat, r_stroke);
        
        int iter = 0;
        for (; iter < _maxIter; iter++)
        {
            // 2c
            // Ap:
            matrix.Mul(p, Ap);
            // alpha:
            var pAp = Dot(p, Ap);
            var alpha = rr0 / pAp;

            // 6d
            Axpy(alpha, p, x);

            // 6e
            Axpy(-alpha, Ap, r_hat);

            // 6f
            r_hat.CopyTo(r_stroke);
            
            matrix.InvLMul(r_stroke);
            Vmul(r_stroke, matrix.Di);
            matrix.InvUMul(r_stroke);

            // 6g
            var rr1 = Dot(r_hat, r_stroke);
            var b_hat = rr1 / rr0;

            // 6h
            Scale(b_hat, p);
            Axpy(1, r_stroke, p);

            rr0 = rr1;

            r_hat.CopyTo(di_inv);
            // matrix.LMul(r_hat, di_inv);
            // matrix.LMul(r_hat, di_inv);
            var rr = Dot(r_hat, r_hat);
            var bb = Dot(b, b);
            if (rr / bb < _eps)
            {
                break;
            }
        }

        matrix.Mul(x, z);
        b.CopyTo(r_hat);
        Axpy(-1, z, r_hat);
        // BLAS.axpy(_x.Length, -1, t, r);
        var rr2 = Dot(r_hat, r_hat);

        return (rr2, iter);
    }
}
