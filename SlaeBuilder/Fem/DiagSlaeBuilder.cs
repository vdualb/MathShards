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

using System.Diagnostics;

namespace MathShards.SlaeBuilder.Fem;

using Dim2 = FiniteElements.Rectangle.Lagrange.BiLinear;
using Dim1 = FiniteElements.Line.Lagrange.Linear;
using MathShards.TelmaCore;
using MathShards.Mesh.RectMesh;
using MathShards.Matrices.Types;
using MathShards.Matrices;
using MathShards.Fem.Common;

using static MathShards.Quadrature.Gauss;

class PatchVectorHAligned : IPatchVector
{
    public required int[] Dofs { get; init; }
    public double[] Values { get => throw new NotImplementedException(); init => throw new NotImplementedException(); }
}

// class PatchVectorVAligned

interface IPatchVector
{
    int[] Dofs { get; init; }
    Real[] Values { get; init; }
}

interface IPatchMatrix
{
    int[] Dofs { get; init; }
    Real[,] Values { get; init; } 
}

public class DiagSlaeBuilder<Tc> : IFemSlaeBuilder
where Tc : CoordSystem.Dim2.ICoordSystem
{
    Diag9Matrix _matrix;
    Real[] _b = [];

    readonly FemRectMesh _mesh;
    public FemRectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly TaskFuncs _funcs;

    public DiagSlaeBuilder(FemRectMesh mesh, TaskFuncs funcs)
    {
        _mesh = mesh;
        _matrix = new Diag9Matrix();
        _funcs = funcs;
    }

    public static IFemSlaeBuilder Construct(FemRectMesh mesh, TaskFuncs funcs)
        => new DiagSlaeBuilder<Tc>(mesh, funcs);

    public (IMatrix, Real[]) Build()
    {
        Trace.WriteLine($"Diag Builder: {GlobalMatrixImpl}");

        Trace.Indent();
        var sw = Stopwatch.StartNew();
        GlobalMatrixInit();
        Trace.WriteLine($"Init: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        GlobalMatrixBuild();
        Trace.WriteLine($"Build: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        BoundaryConditionsApply();
        Trace.WriteLine($"Conds: {sw.ElapsedMilliseconds}ms");
        Trace.Unindent();

        // HACK: Добавка в правую часть от точечного источника
        // var idx = _mesh.GetDofAtInitNode(0, 2);
        // _b[idx] += 1 / (2 * Math.PI);

        return (_matrix, _b);
    }
    
    void GlobalMatrixPortraitCompose()
    {
        // в случае диагональной матрицы это просто Gap
        _matrix.Gap = _mesh.X.Length - 2;
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _matrix.Di = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _matrix.Ld0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _matrix.Rd0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _b = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
    }

    void BoundaryConditionType1Apply(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        void SharedBody(int m, Real b)
        {
            _b[m] = b;
            _matrix.Di[m] = 1;

            int t;

            /* Гауссово исключение столбца */
            t = m-3-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.Rd3[t];
            t = m-2-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.Rd2[t];
            t = m-1-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.Rd1[t];
            t = m-1;
            if (t >= 0) _b[t] -= b * _matrix.Rd0[t];
            
            t = m+1;
            if (t < _matrix.Size) _b[t] -= b * _matrix.Ld0[m];
            t = m+1+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.Ld1[m];
            t = m+2+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.Ld2[m];
            t = m+3+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.Ld3[m];

            /* Обнуление строки и столбца */
            _matrix.Rd3[m] = 0;
            _matrix.Rd2[m] = 0;
            _matrix.Rd1[m] = 0;
            _matrix.Rd0[m] = 0;
            _matrix.Ld3[m] = 0;
            _matrix.Ld2[m] = 0;
            _matrix.Ld1[m] = 0;
            _matrix.Ld0[m] = 0;

            t = m - 3 - _matrix.Gap;
            if (t >= 0)
            {
                _matrix.Ld3[t] = 0;
                _matrix.Rd3[t] = 0;
            }
            t = m - 2 - _matrix.Gap;
            if (t >= 0) 
            {
                _matrix.Ld2[t] = 0;
                _matrix.Rd2[t] = 0;
            }
            t = m - 1 - _matrix.Gap;
            if (t >= 0) 
            {
                _matrix.Ld1[t] = 0;
                _matrix.Rd1[t] = 0;
            }
            t = m - 1;
            if (t >= 0) 
            {
                _matrix.Ld0[t] = 0;
                _matrix.Rd0[t] = 0;
            }
        }
        
        if (x1 == x2)
        {
            for (int yi = y1; yi <= y2; yi++)
            {
                var m = yi * _mesh.X.Length + x1;
                var b = _funcs.Ug(num, _mesh.X[x1], _mesh.Y[yi]);
                SharedBody(m, b);
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi <= x2; xi++)
            {
                var m = y1 * _mesh.X.Length + xi;
                var b = _funcs.Ug(num, _mesh.X[xi], _mesh.Y[y1]);
                SharedBody(m, b);
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    // HACK: две копипасты impl
    void Cond2Impl(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int xi1 = _mesh.XAfterGridInit(bc.X1);
        int xi2 = _mesh.XAfterGridInit(bc.X2);
        int yi1 = _mesh.YAfterGridInit(bc.Y1);
        int yi2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */
        
        var num = bc.Num;

        if (xi1 == xi2)
        {
            // вертикальная граница
            for (int yi = yi1; yi < yi2; yi++)
            {
                Real x = _mesh.X[xi1];
                Real y0 = _mesh.Y[yi];
                Real y1 = _mesh.Y[yi + 1];
                var hy = y1 - y0;

                var n = new int[2];
                n[0] = yi * _mesh.X.Length + xi1;
                n[1] = (yi + 1) * _mesh.X.Length + xi2;
                
                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    _b[n[i]] += Integrate1DOrder5(y0, y1, 
                        y =>
                        {
                            // в координатах шаблонного базиса - [0;1]
                            var y01 = (y - y0) / hy;
                            return _funcs.Theta(num, x, y)
                                * Dim1.Basis[i](y01)
                                * Tc.Jacobian(x, y);
                        }
                    );
                }
            }
        }
        else if (yi1 == yi2)
        {
            // горизонтальная граница
            for (int xi = xi1; xi < xi2; xi++)
            {
                Real y = _mesh.Y[yi1];
                Real x0 = _mesh.X[xi];
                Real x1 = _mesh.X[xi + 1];
                var hx = x1 - x0;

                int[] n = [
                    yi1 * _mesh.X.Length + xi,
                    yi2 * _mesh.X.Length + xi + 1    
                ];

                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    _b[n[i]] += Integrate1DOrder5(x0, x1, 
                        x =>
                        {
                            // в координатах шаблонного базиса - [0;1]
                            var x01 = (x - x0) / hx;
                            return _funcs.Theta(num, x, y)
                                * Dim1.Basis[i](x01)
                                * Tc.Jacobian(x, y);
                        }
                    );
                    
                }
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }
    
    void Cond2ImplOld(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int xi1 = _mesh.XAfterGridInit(bc.X1);
        int xi2 = _mesh.XAfterGridInit(bc.X2);
        int yi1 = _mesh.YAfterGridInit(bc.Y1);
        int yi2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */
        
        var num = bc.Num;

        if (xi1 == xi2)
        {
            // вертикальная граница
            for (int yi = yi1; yi < yi2; yi++)
            {
                var h = _mesh.Y[yi + 1] - _mesh.Y[yi];

                Real k1 = _funcs.Theta(num, _mesh.X[xi1], _mesh.Y[yi]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[xi2], _mesh.Y[yi + 1]);
                int n1 = yi * _mesh.X.Length + xi1;
                int n2 = (yi + 1) * _mesh.X.Length + xi2;
                
                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else if (yi1 == yi2)
        {
            // горизонтальная граница
            for (int xi = xi1; xi < xi2; xi++)
            {
                var h = _mesh.X[xi + 1] - _mesh.X[xi];

                Real k1 = _funcs.Theta(num, _mesh.X[xi], _mesh.Y[yi1]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[xi + 1], _mesh.Y[yi2]);
                int n1 = yi1 * _mesh.X.Length + xi;
                int n2 = yi2 * _mesh.X.Length + xi + 1;

                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }
    
    void BoundaryConditionType2Apply(BoundaryCondition bc)
    {
        Cond2Impl(bc);
    }

    void Cond3ImplOld(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2]; // 'hat B'
        var localA = new Real[2, 2]; // 'hat A'
        Real h;
        if (x1 == x2)
        {
            for (int yi = y1; yi < y2; yi++)
            {
                h = _mesh.Y[yi + 1] - _mesh.Y[yi];
                localA[0, 0] = localA[1, 1] = _funcs.Beta(num) * h / 3;
                localA[0, 1] = localA[1, 0] = _funcs.Beta(num) * h / 6;

                Real k1 = _funcs.uBeta(num, _mesh.X[x1], _mesh.Y[yi]);
                Real k2 = _funcs.uBeta(num, _mesh.X[x2], _mesh.Y[yi + 1]);
                localB[0] = h * _funcs.Beta(num) * (2 * k1 + k2) / 6;
                localB[1] = h * _funcs.Beta(num) * (k1 + 2 * k2) / 6;

                var m = new int[2];
                m[0] = yi * _mesh.X.Length + x1;
                m[1] = (yi + 1) * _mesh.X.Length + x2;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.Rd2[m[0]] += localA[0, 1];
                _matrix.Ld2[m[0]] += localA[1, 0];
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi < x2; xi++)
            {
                h = _mesh.X[xi + 1] - _mesh.X[xi];
                localA[0, 0] = localA[1, 1] = _funcs.Beta(num) * h / 3;
                localA[0, 1] = localA[1, 0] = _funcs.Beta(num) * h / 6;

                Real k1 = _funcs.uBeta(num, _mesh.X[xi], _mesh.Y[y1]);
                Real k2 = _funcs.uBeta(num, _mesh.X[xi + 1], _mesh.Y[y2]);
                localB[0] = h * _funcs.Beta(num) * (2 * k1 + k2) / 6;
                localB[1] = h * _funcs.Beta(num) * (k1 + 2 * k2) / 6;

                var m = new int[2];
                m[0] = y1 * _mesh.X.Length + xi;
                m[1] = y2 * _mesh.X.Length + xi + 1;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.Rd0[m[0]] += localA[0, 1];
                _matrix.Ld0[m[0]] += localA[1, 0];
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }
    
    void Cond3Impl(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int xi1 = _mesh.XAfterGridInit(bc.X1);
        int xi2 = _mesh.XAfterGridInit(bc.X2);
        int yi1 = _mesh.YAfterGridInit(bc.Y1);
        int yi2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2]; // 'hat B'
        var localA = new Real[2, 2]; // 'hat A'
        if (xi1 == xi2)
        {
            for (int yi = yi1; yi < yi2; yi++)
            {
                var y0 = _mesh.Y[yi];
                var y1 = _mesh.Y[yi + 1];
                var x = _mesh.X[xi1];

                var hy = _mesh.Y[yi + 1] - _mesh.Y[yi];

                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    for (int j = 0; j < Dim1.Basis.Length; j++)
                    {
                        localA[i, j] = Integrate1DOrder5(y0, y1,
                            y =>
                            {
                                // в координатах шаблонного базиса - [0;1]
                                var y01 = (y - y0) / hy;
                                return _funcs.Beta(num)
                                    * Dim1.Basis[i](y01)
                                    * Dim1.Basis[j](y01)
                                    * Tc.Jacobian(x, y);
                            }
                        );
                    }
                }

                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    localB[i] = Integrate1DOrder5(y0, y1, 
                        y =>
                        {
                            // в координатах шаблонного базиса - [0;1]
                            var y01 = (y - y0) / hy;
                            return _funcs.Beta(num)
                                * _funcs.uBeta(num, x, y)
                                * Dim1.Basis[i](y01)
                                * Tc.Jacobian(x, y);
                        }
                    );
                }
                
                var m = new int[2];
                m[0] = yi * _mesh.X.Length + xi1;
                m[1] = (yi + 1) * _mesh.X.Length + xi2;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.Rd2[m[0]] += localA[0, 1];
                _matrix.Ld2[m[0]] += localA[1, 0];
            }
        }
        else if (yi1 == yi2)
        {
            // по горизонтали
            // TODO: даёт неверное решение
            throw new NotImplementedException("Третьи условия сломаны");
            for (int xi = xi1; xi < xi2; xi++)
            {
                var x0 = _mesh.X[xi];
                var x1 = _mesh.X[xi + 1];
                var y = _mesh.Y[yi1];

                var hx = _mesh.X[xi + 1] - _mesh.X[xi];

                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    for (int j = 0; j < Dim1.Basis.Length; j++)
                    {
                        localA[i, j] = Integrate1DOrder5(x0, x1,
                            x =>
                            {
                                // в координатах шаблонного базиса - [0;1]
                                var x01 = (x - x0) / hx;
                                return _funcs.Beta(num)
                                    * Dim1.Basis[i](x01)
                                    * Dim1.Basis[j](x01)
                                    * x;
                            }
                        );
                    }
                }

                for (int i = 0; i < Dim1.Basis.Length; i++)
                {
                    localB[i] = Integrate1DOrder5(x0, x1,
                        x =>
                        {
                            // в координатах шаблонного базиса - [0;1]
                            var x01 = (x - x0) / hx;
                            return _funcs.Beta(num)
                                * _funcs.uBeta(num, x, y)
                                * Dim1.Basis[i](x01)
                                * x;
                        }
                    );
                }

                var m = new int[2];
                m[0] = yi1 * _mesh.X.Length + xi;
                m[1] = yi2 * _mesh.X.Length + xi + 1;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.Rd2[m[0]] += localA[0, 1];
                _matrix.Ld2[m[0]] += localA[1, 0];
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }
    
    void BoundaryConditionType3Apply(BoundaryCondition bc)
    {
        Cond3Impl(bc);
    }

    void BoundaryConditionsApply()
    {
        var bc_type1 = new List<BoundaryCondition>();

        foreach (var bc in _mesh.BoundaryConditions)
        {
            var type = bc.Type;

            switch (type)
            {
                case 1:
                    /* К.у. первого рода будут применены последними */
                    bc_type1.Add(bc);
                    break;
                case 2:
                    BoundaryConditionType2Apply(bc);
                    break;
                case 3:
                    BoundaryConditionType3Apply(bc);
                    break;

                default:
                    throw new Exception("Странный тип краевого условия");
            }
        }

        foreach (var b1 in bc_type1)
        {
            BoundaryConditionType1Apply(b1);
        }
    }

    void GlobalMatrixBuild()
    {
        switch (GlobalMatrixImpl)
        {
            case GlobalMatrixImplType.Host:
                GlobalMatrixBuildImplHost();
                break;
            default:
                throw new InvalidOperationException();
        }
    }

    void GlobalMatrixBuildImplHost()
    {
        for (int yi = 0; yi < _mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
            {
                var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);
                if (subDom == null) continue;
                
                PairF64 p0 = new(_mesh.X[xi], _mesh.Y[yi]);
                PairF64 p1 = new(_mesh.X[xi + 1], _mesh.Y[yi + 1]);
                
                var local = Dim2.ComputeLocal<Tc>(_funcs, p0, p1, subDom.Value);

                int a = yi * _mesh.X.Length + xi;

                _matrix.Di[a] += local[0, 0];
                _matrix.Ld0[a] += local[0, 1];
                _matrix.Rd0[a] += local[0, 1];
                _matrix.Ld2[a] += local[0, 2];
                _matrix.Rd2[a] += local[0, 2];
                _matrix.Ld3[a] += local[0, 3];
                _matrix.Rd3[a] += local[0, 3];

                _matrix.Di[a + 1] += local[1, 1];
                _matrix.Ld1[a + 1] += local[1, 2];
                _matrix.Rd1[a + 1] += local[1, 2];
                _matrix.Ld2[a + 1] += local[1, 3];
                _matrix.Rd2[a + 1] += local[1, 3];

                var a2 = a + _mesh.X.Length;
                _matrix.Di[a2] += local[2, 2];
                _matrix.Ld0[a2] += local[2, 3];
                _matrix.Rd0[a2] += local[2, 3];

                _matrix.Di[a2 + 1] += local[3, 3];

                var localB = Dim2.ComputeLocalB<Tc>(_funcs, p0, p1, subDom.Value);

                _b[a] += localB[0];
                _b[a + 1] += localB[1];
                _b[a2] += localB[2];
                _b[a2 + 1] += localB[3];
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _matrix.Di.Length; i++)
        {
            if (_matrix.Di[i] == 0)
            {
                _matrix.Di[i] = 1;
            }
        }
    }
}
