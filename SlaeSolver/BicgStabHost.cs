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

namespace MathShards.SlaeSolver;
using static MathShards.SlaeSolver.Shared;

// в какой-то момент этот класс превратился в BiCGStabMkl
public class BicgStabHost : ISlaeSolver
{
    int _maxIter;
    Real _eps;

    int _n = 0; // размерность СЛАУ
    Real[] r = [];
    Real[] di_inv = [];
    Real[] y = [];
    Real[] z = [];
    Real[] ks = [];
    Real[] kt = [];
    Real[] r_hat = [];
    Real[] p = [];
    Real[] nu = [];
    Real[] h = [];
    Real[] s = [];
    Real[] t = [];

    public BicgStabHost(int maxIter, Real eps)
    {
        _maxIter = maxIter;
        _eps = eps;
    }

    public static ISlaeSolver Construct(int maxIter, Real eps)
        => new BicgStabHost(maxIter, eps);

    // Выделить память для временных массивов
    // n - длина каждого массива
    public void AllocateTemps(int n)
    {
        if (n != _n)
        {
            _n = n;

            r = new Real[_n];
            r_hat = new Real[_n];
            p = new Real[_n];
            nu = new Real[_n];
            h = new Real[_n];
            s = new Real[_n];
            t = new Real[_n];
            di_inv = new Real[_n];
            y = new Real[_n];
            z = new Real[_n];
            ks = new Real[_n];
            kt = new Real[_n];
        }
    }

    public (Real discrep, int iter) Solve<T>(T matrix, Span<Real> b, Span<Real> x)
    where T: Matrices.Types.IMatrix
    {
        AllocateTemps(x.Length);

        var _b = b;

        var r = this.r.AsSpan();
        var r_hat = this.r_hat.AsSpan();
        var p = this.p.AsSpan();
        var nu = this.nu.AsSpan();
        var h = this.h.AsSpan();
        var s = this.s.AsSpan();
        var t = this.t.AsSpan();
        var di_inv = this.di_inv.AsSpan();
        var y = this.y.AsSpan();
        var z = this.z.AsSpan();
        var ks = this.ks.AsSpan();
        var kt = this.kt.AsSpan();


        // precond
        matrix.Di.CopyTo(di_inv);
        Rsqrt(di_inv);
        // 1.
        matrix.Mul(x, t);
        _b.CopyTo(r);
        Axpy(-1, t, r);
        // 2.
        r.CopyTo(r_hat);
        // 3.
        Real pp = Dot(r, r); // r_hat * r
        // 4.
        r.CopyTo(p);

        int iter = 0;
        Real rr;
        for (; iter < _maxIter; iter++)
        {
            // 1.
            p.CopyTo(y);
            Vmul(y, di_inv);
            Vmul(y, di_inv);

            // 2.
            matrix.Mul(y, nu);

            // 3.
            Real rnu = Dot(r_hat, nu);
            Real alpha = pp / rnu;

            // 4.
            x.CopyTo(h);
            Axpy(alpha, y, h);
            // BLAS.axpy(_n, alpha, y, h);

            // 5.
            r.CopyTo(s);
            Axpy(-alpha, nu, s);
            // BLAS.axpy(_n, -alpha, nu, s);

            // 6.
            Real ss = Dot(s, s);
            if (ss < _eps)
            {
                h.CopyTo(x);
                // _x.Dispose();
                // _x = h;
                break;
            }

            // 7.
            s.CopyTo(ks);
            Vmul(ks, di_inv);
            ks.CopyTo(z);
            Vmul(z, di_inv);

            // 8.
            matrix.Mul(z, t);

            // 9.
            t.CopyTo(kt);
            Vmul(kt, di_inv);

            Real ts = Dot(ks, kt);
            Real tt = Dot(kt, kt);
            Real w = ts / tt;

            // 10.
            h.CopyTo(x);
            Axpy(w, z, x);
            // BLAS.axpy(_n, w, z, _x);

            // 11.
            s.CopyTo(r);
            Axpy(-w, t, r);
            // BLAS.axpy(_n, -w, t, r);

            // 12.
            rr = Dot(r, r);
            if (rr < _eps)
            {
                break;
            }

            // 13-14
            Real pp1 = Dot(r, r_hat);
            Real beta = (pp1 / pp) * (alpha / w);

            // 15.
            Axpy(-w, nu, p);
            // BLAS.axpy(_n, -w, nu, p);
            Scale(beta, p);
            // BLAS.scal(_n, beta, p);
            // BLAS.axpy(_n, 1, r, p);
            Axpy(1, r, p);

            pp = pp1;
        }

        matrix.Mul(x, t);
        _b.CopyTo(r);
        Axpy(-1, t, r);
        // BLAS.axpy(_x.Length, -1, t, r);
        rr = Dot(r, r);

        return (rr, iter);
    }
}
