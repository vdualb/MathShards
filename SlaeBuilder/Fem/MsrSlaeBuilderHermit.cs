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

using System.Diagnostics;

namespace MathShards.SlaeBuilder.Fem;

using MathShards.TelmaCore;
using MathShards.Mesh.RectMesh;
using MathShards.Matrices.Types;
using MathShards.Matrices;
using MathShards.Fem.Common;
using Dim2 = FiniteElements.Rectangle.Hermit.Cubic;

class DiagSlaeBuilderHermit<Tc> : IFemSlaeBuilder
where Tc : CoordSystem.Dim2.ICoordSystem
{
    MsrMatrix _matrix;
    Real[] _b = [];

    readonly FemRectMesh _mesh;
    public FemRectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly ITaskFuncs _funcs;

    public DiagSlaeBuilderHermit(FemRectMesh mesh, ITaskFuncs funcs)
    {
        _mesh = mesh;
        _matrix = new MsrMatrix();
        _funcs = funcs;
    }

    public static IFemSlaeBuilder Construct(FemRectMesh mesh, ITaskFuncs funcs)
        => new DiagSlaeBuilderHermit<Tc>(mesh, funcs);

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
        Trace.WriteLine($"Conds: {sw.ElapsedMilliseconds}");
        Trace.Unindent();

        return (_matrix, _b);
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        int n = _matrix.Ia.Length - 1;

        _matrix.Di = Enumerable.Repeat((Real)0, n).ToArray();
        _matrix.Elems = Enumerable.Repeat((Real)0, _matrix.Ja.Length).ToArray();

        _b = Enumerable.Repeat((Real)0, n).ToArray();
    }

    void GlobalMatrixPortraitCompose()
    {
        int dofsPerNode = 4;
        // количество узлов, лежащих на элементе
        int numberOfUnknowns(int i) => dofsPerNode * 4;
        // (номер элемента, локальный номер узла) -> глобальный номер узла
        int idxOfUnknown(int ielem, int j)
        {
            int x0 = ielem % (_mesh.X.Length - 1);
            int y0 = ielem / (_mesh.X.Length - 1);

            if (j < 8)
            {
                return y0*_mesh.X.Length*dofsPerNode + x0*dofsPerNode + j;
            }
            else if (j < 16)
            {
                return (y0 + 1)*_mesh.X.Length*dofsPerNode + x0*dofsPerNode + (j - dofsPerNode*2);
            }
            else
            {
                throw new ArgumentException("Странная координата конечного элемента");
            }
        }

        HashSet<int>[] list = new HashSet<int>[Mesh.nodesCount*dofsPerNode];
        for (int i = 0; i < list.Length; i++)
        {
            list[i] = [];
        }

        /* цикл по всем конечным элементам */
        for (int ielem = 0; ielem < Mesh.feCount; ielem++)
        {
            /* цикл по всем узлам данного к.э. */
            for (int idx0 = 0; idx0 < numberOfUnknowns(ielem); idx0++)
            {
                /* цикл по узлам, соседним с idx0 */
                for (int idx1 = 0; idx1 < numberOfUnknowns(ielem); idx1++)
                {
                    if (idx0 == idx1) continue;
                    /* нахождение глобальных номеров локальных узлов */
                    int k1 = idxOfUnknown(ielem, idx0);
                    int k2 = idxOfUnknown(ielem, idx1);
                    /* */
                    /* заносим в list[k2] номера всех узлов, "соседних" с ним.
                        По сути здесь k2 - номер строки, а k1 - номер ненулевого
                        элемента в глобальной матрице*/
                    list[k1].Add(k2);
                }
            }
        }

        _matrix.Ia = new int[list.Length + 1];
        _matrix.Ia[0] = 0;
        /* формирование массивов ig jg по списку list */
        for (int i = 1; i < _matrix.Ia.Length; i++)
        {
            _matrix.Ia[i] = _matrix.Ia[i - 1] + list[i - 1].Count;
        }
        _matrix.Ja = new int[_matrix.Ia[list.Length]];
        for (var i = 0; i < list.Length; i++)
        {
            var row = list[i].Order().ToArray();
            for (int j = _matrix.Ia[i]; j < _matrix.Ia[i + 1]; j++)
            {
                _matrix.Ja[j] = row[j - _matrix.Ia[i]];
            }
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
        int dofsPerNode = 4;
        for (int yi = 0; yi < _mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
            {
                var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);
                if (subDom == null) continue;

                PairReal p0 = new(_mesh.X[xi], _mesh.Y[yi]);
                PairReal p1 = new(_mesh.X[xi + 1], _mesh.Y[yi + 1]);

                var local = Dim2.ComputeLocal<Tc>(_funcs, p0, p1, subDom.Value);
                var localB = Dim2.ComputeLocalB<Tc>(_funcs, p0, p1, subDom.Value);
                
                var dofs = new int[16];
                dofs[0] = yi*_mesh.X.Length*dofsPerNode + xi*dofsPerNode;
                for (int i = 1; i < 8; i++)
                {
                    dofs[i] = dofs[i-1]+1;
                }
                dofs[8] = (yi + 1)*_mesh.X.Length*dofsPerNode + xi*dofsPerNode - dofsPerNode*2;
                for (int i = 9; i < 16; i++)
                {
                    dofs[i] = dofs[i-1]+1;
                }
                
                for (int i = 0; i < 16; i++)
                {
                    int a = _matrix.Ia[dofs[i]];
                    for (int j = 0; j < 16; j++)
                    {
                        if (i == j)
                        {
                            _matrix.Di[dofs[i]] = local[i, j];
                        } else {
                            a = Shared.LFind(_matrix.Ja, dofs[j], a);
                        }
                    }
                    _b[dofs[i]] += localB[i];
                }
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
    
    void BoundaryConditionType1Apply(BoundaryCondition bc)
    {
        throw new NotImplementedException();
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        if (x1 == x2)
        {
            for (int yi = y1; yi <= y2; yi++)
            {
                var m = yi * _mesh.X.Length + x1;
                _b[m] = _funcs.Ug(num, _mesh.X[x1], _mesh.Y[yi]);
                _matrix.Di[m] = 1;

                /* Обнуление строки */
                int ig0 = _matrix.Ia[m];
                int ig1 = _matrix.Ia[m + 1];
                for (int j = ig0; j < ig1; j++)
                {
                    _matrix.Elems[j] = 0;
                }
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi <= x2; xi++)
            {
                var m = y1 * _mesh.X.Length + xi;
                _b[m] = _funcs.Ug(num, _mesh.X[xi], _mesh.Y[y1]);
                _matrix.Di[m] = 1;

                /* Обнуление строки */
                int ig0 = _matrix.Ia[m];
                int ig1 = _matrix.Ia[m + 1];
                for (int j = ig0; j < ig1; j++)
                {
                    _matrix.Elems[j] = 0;
                }
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    void BoundaryConditionType2Apply(BoundaryCondition bc)
    {
        throw new NotImplementedException();
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        if (x1 == x2)
        {
            for (int yi = y1; yi < y2; yi++)
            {
                var h = _mesh.Y[yi + 1] - _mesh.Y[yi];

                Real k1 = _funcs.Theta(num, _mesh.X[x1], _mesh.Y[yi]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[x2], _mesh.Y[yi + 1]);
                int n1 = yi * _mesh.X.Length + x1;
                int n2 = (yi + 1) * _mesh.X.Length + x2;

                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi < x2; xi++)
            {
                var h = _mesh.X[xi + 1] - _mesh.X[xi];

                Real k1 = _funcs.Theta(num, _mesh.X[xi], _mesh.Y[y1]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[xi + 1], _mesh.Y[y2]);
                int n1 = y1 * _mesh.X.Length + xi;
                int n2 = y2 * _mesh.X.Length + xi + 1;

                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    void BoundaryConditionType3Apply(BoundaryCondition bc)
    {
        throw new NotImplementedException();
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

                var res = Shared.LFind(_matrix.Ja, what: m[1], start: _matrix.Ia[m[0]]);
                _matrix.Elems[res] += localA[0, 1];
                res = Shared.LFind(_matrix.Ja, what: m[0], start: _matrix.Ia[m[1]]);
                _matrix.Elems[res] += localA[1, 0];
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

                var res = Shared.LFind(_matrix.Ja, what: m[1], start: _matrix.Ia[m[0]]);
                _matrix.Elems[res] += localA[0, 1];
                res = Shared.LFind(_matrix.Ja, what: m[0], start: _matrix.Ia[m[1]]);
                _matrix.Elems[res] += localA[1, 0];
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }
}
