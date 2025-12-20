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

namespace MathShards.SlaeBuilder.Spline;

using Matrices;
using Dim1 = FiniteElements.Line.Hermit.Cubic;
using static Quadrature.Gauss;
using static SlaeBuilder.Fem.Shared;
using Mesh.RectMesh;
using Matrices.Types;

public class MsrSlaeBuilderHermit1D : ISplineSlaeBuilder1D
{
    MsrMatrix _matrix;
    Real[] _b = [];
    
    // сетка, на которой лежат значения _inValues
    readonly LineMesh _inMesh;
    readonly Real[] _inValues;
    
    // сетка, на которой строится сплайн
    readonly LineMesh _splineMesh;
    readonly SplineParams _splineParams;

    public MsrSlaeBuilderHermit1D(
        LineMesh inMesh, Real[] inValues,
        LineMesh splineMesh, SplineParams splineParams
    ) {
        _inMesh = inMesh;
        _inValues = inValues;
        _splineMesh = splineMesh;
        _splineParams = splineParams;
        
        _matrix = new MsrMatrix();
    }
    
    public static ISplineSlaeBuilder1D Construct(
        LineMesh inMesh, double[] inValues,
        LineMesh mesh, SplineParams splineParams
    ) {
        return new MsrSlaeBuilderHermit1D(
            inMesh, inValues,
            mesh, splineParams
        );
    }

    public (IMatrix matrix, double[] right) Build()
    {
        Trace.Indent();
        var sw = Stopwatch.StartNew();
        GlobalMatrixInit();
        Trace.WriteLine($"Init: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        GlobalMatrixBuild();
        Trace.WriteLine($"Build: {sw.ElapsedMilliseconds}ms");

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
        int dofsPerNode = 2;
        int nodesPerElement = 2;
        // количество узлов, лежащих на элементе
        int numberOfUnknowns(int i) => dofsPerNode * nodesPerElement;
        // (номер элемента, локальный номер узла) -> глобальный номер узла
        int idxOfUnknown(int ielem, int j)
        {
            int start = ielem * nodesPerElement;
            
            return start + j;
        }

        HashSet<int>[] list = new HashSet<int>[_splineMesh.nodesCount*dofsPerNode];
        for (int i = 0; i < list.Length; i++)
        {
            list[i] = [];
        }

        var feCount = _splineMesh.feCount;
        /* цикл по всем конечным элементам */
        for (int ielem = 0; ielem < feCount; ielem++)
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

    Real[,] ComputeLocal(Real p0, Real p1, Real srcp)
    {
        // TODO: способ задания w извне
        var w = _splineParams.W;
        var side = 4;
        var результат = new Real[side, side];
        
        for (int i = 0; i < side; i++)
        {
            for (int j = 0; j < side; j++)
            {
                var val = Dim1.BasisConverted(i, p0, p1, srcp)
                    * Dim1.BasisConverted(j, p0, p1, srcp);
                результат[i, j] = w * val;
            }
        }

        return результат;
    }
    
    Real[,] ComputeLocalRegDiff1(Real p0, Real p1)
    {
        var alpha = _splineParams.Alpha;
        var side = 4;
        var результат = new Real[side, side];
        
        for (int i = 0; i < side; i++)
        {
            for (int j = 0; j < side; j++)
            {
                var func = (Real point) =>
                {
                    return Dim1.BasisGradConverted(i, p0, p1, point)
                    * Dim1.BasisGradConverted(j, p0, p1, point);
                };

                результат[i, j] = alpha * Integrate1DOrder5(p0, p1, func);
            }
        }

        return результат;
    }
    
    Real[,] ComputeLocalRegDiff2(Real p0, Real p1)
    {
        var beta = _splineParams.Beta;
        var side = 4;
        var результат = new Real[side, side];
        
        for (int i = 0; i < side; i++)
        {
            for (int j = 0; j < side; j++)
            {
                var func = (Real point) =>
                {
                    return Dim1.BasisGradGradConverted(i, p0, p1, point)
                    * Dim1.BasisGradGradConverted(j, p0, p1, point);
                };

                результат[i, j] = beta * Integrate1DOrder5(p0, p1, func);
            }
        }

        return результат;
    }
    
    Real[] ComputeLocalB(Real p0, Real p1, Real srcp, int dof)
    {
        var w = _splineParams.W;
        var side = 4;
        var результат = new Real[side];
        
        for (int i = 0; i < side; i++)
        {
            результат[i] = w * Dim1.BasisConverted(i, p0, p1, srcp) * _inValues[dof];
        }

        return результат;
    }
    
    void AddLocal(int[] dofs, Real[,] mat) {
        for (int i = 0; i < 4; i++)
        {
            int a = _matrix.Ia[dofs[i]];
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                {
                    _matrix.Di[dofs[i]] += mat[i, j];
                } else {
                    a = LFind(_matrix.Ja, dofs[j], a);
                    _matrix.Elems[a] += mat[i, j];
                }
            }
        }
    }
    
    void AddLocal(int[] dofs, Real[] locVec) {
        for (int i = 0; i < 4; i++)
        {
            _b[dofs[i]] += locVec[i];
        }
    }
    
    // TODO: починить для случая с разбиениями
    void GlobalMatrixBuild()
    {
        int dofsPerNode = 2;
        
        // обход по всем сплайнам
        int srci0 = 0;
        for (int spli = 0; spli < _splineMesh.X.Length - 1; spli++)
        {
            Real p0 = _splineMesh.X[spli];
            Real p1 = _splineMesh.X[spli + 1];

            int[] dofs = new int[4];
            
            for (int i = 0; i < 4; i++)
            {
                dofs[i] = spli*dofsPerNode + i;
            }

            var reg1 = ComputeLocalRegDiff1(p0, p1);
            var reg2 = ComputeLocalRegDiff2(p0, p1);
            AddLocal(dofs, reg1);
            AddLocal(dofs, reg2);
            
            // обход по конечным элементам внутри сплайна
            int srci = srci0;
            for (
                srci = srci0;
                srci < _inMesh.X.Length && _inMesh.X[srci] <= p1;
                srci++
            ) {
                Real srcp = _inMesh.X[srci];
                
                var local = ComputeLocal(p0, p1, srcp);
                // TODO: verify that it works
                //
                var localB = ComputeLocalB(p0, p1, srcp, srci);
                
                AddLocal(dofs, local);
                AddLocal(dofs, localB);
            }
            srci0 = srci;
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
