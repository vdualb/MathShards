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

namespace MathShards.SlaeBuilder.Spline;

using TelmaCore;
using Mesh.RectMesh;
using Matrices.Types;
using Matrices;
using static Quadrature.Gauss;
using SlaeBuilder.Fem;

using Dim2 = FiniteElements.Rectangle.Hermit.Cubic;


public class MsrSlaeBuilderHermit<Tc> : ISplineSlaeBuilder
where Tc : CoordSystem.Dim2.ICoordSystem
{
    MsrMatrix _matrix;
    Real[] _b = [];
    
    // сетка, на которой лежат значения _inValues
    readonly RectMesh _inMesh;
    readonly Real[] _inValues;
    
    // сетка, на которой строится сплайн
    readonly RectMesh _splineMesh;

    public MsrSlaeBuilderHermit(
        RectMesh inMesh, Real[] inValues,
        RectMesh mesh
    ) {
        _inMesh = inMesh;
        _inValues = inValues;
        _splineMesh = mesh;
        _matrix = new MsrMatrix();
    }

    public static ISplineSlaeBuilder Construct(
        RectMesh inMesh, Real[] inValues,
        RectMesh mesh
    ) {
        return new MsrSlaeBuilderHermit<Tc>(inMesh, inValues, mesh);
    }

    public (IMatrix matrix, Real[] right) Build()
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
        int dofsPerNode = 4;
        // количество узлов, лежащих на элементе
        int numberOfUnknowns(int i) => dofsPerNode * 4;
        // (номер элемента, локальный номер узла) -> глобальный номер узла
        int idxOfUnknown(int ielem, int j)
        {
            int x0 = ielem % (_splineMesh.X.Length - 1);
            int y0 = ielem / (_splineMesh.X.Length - 1);

            if (j < 8)
            {
                return y0*_splineMesh.X.Length*dofsPerNode + x0*dofsPerNode + j;
            }
            else if (j < 16)
            {
                return (y0 + 1)*_splineMesh.X.Length*dofsPerNode + x0*dofsPerNode + (j - dofsPerNode*2);
            }
            else
            {
                throw new ArgumentException("Странная координата конечного элемента");
            }
        }

        HashSet<int>[] list = new HashSet<int>[_splineMesh.nodesCount*dofsPerNode];
        for (int i = 0; i < list.Length; i++)
        {
            list[i] = [];
        }

        /* цикл по всем конечным элементам */
        for (int ielem = 0; ielem < _splineMesh.feCount; ielem++)
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

    Real[,] ComputeLocal(PairReal p0, PairReal p1, PairReal srcp)
    {
        // TODO: способ задания w извне
        var w = (Real)0.9;
        var результат = new Real[16, 16];
        
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 16; j++)
            {
                var val = Dim2.BasisConverted(i, p0, p1, srcp)
                    * Dim2.BasisConverted(j, p0, p1, srcp);
                результат[i, j] = w * val;
            }
        }

        return результат;
    }
    
    Real[] ComputeLocalB(PairReal p0, PairReal p1, PairReal srcp, int dof)
    {
        // TODO: способ задания w извне
        var w = (Real)0.9;
        var результат = new Real[16];
        
        for (int i = 0; i < 16; i++)
        {
            результат[i] = w * Dim2.BasisConverted(i, p0, p1, srcp) * _inValues[dof];
        }

        return результат;
    }
    
    void GlobalMatrixBuild()
    {
        int dofsPerNode = 4;
        
        // обход по всем сплайнам
        int srcyi0 = 0;
        for (int yi = 0; yi < _splineMesh.Y.Length - 1; yi++)
        {
            int srcxi0 = 0;
            int srcyi = srcyi0;
            for (int xi = 0; xi < _splineMesh.X.Length - 1; xi++)
            {
                PairReal p0 = new(_splineMesh.X[xi], _splineMesh.Y[yi]);
                PairReal p1 = new(_splineMesh.X[xi + 1], _splineMesh.Y[yi + 1]);

                var dofs = new int[16];
                {
                    dofs[0] = yi*_splineMesh.X.Length*dofsPerNode + xi*dofsPerNode;
                    for (int i = 1; i < 8; i++)
                    {
                        dofs[i] = dofs[i-1]+1;
                    }
                    dofs[8] = (yi + 1)*_splineMesh.X.Length*dofsPerNode + xi*dofsPerNode;
                    for (int i = 9; i < 16; i++)
                    {
                        dofs[i] = dofs[i-1]+1;
                    }   
                }
                
                // обход по конечным элементам внутри сплайна
                int srcxi = srcxi0;
                {
                    for (
                        srcyi = srcyi0;
                        srcyi < _inMesh.Y.Length && _inMesh.Y[srcyi] <= p1.Y;
                        srcyi++
                    ) {
                        for (
                            srcxi = srcxi0;
                            srcxi < _inMesh.X.Length && _inMesh.X[srcxi] <= p1.X;
                            srcxi++
                        ) {
                            PairReal srcp = new(_inMesh.X[srcxi], _inMesh.Y[srcyi]);
                            var local = ComputeLocal(p0, p1, srcp);
                            int dof = srcyi * _inMesh.X.Length + srcxi;
                            var localB = ComputeLocalB(p0, p1, srcp, dof);
                            
                            for (int i = 0; i < 16; i++)
                            {
                                int a = _matrix.Ia[dofs[i]];
                                for (int j = 0; j < 16; j++)
                                {
                                    if (i == j)
                                    {
                                        _matrix.Di[dofs[i]] += local[i, j];
                                    } else {
                                        a = Shared.LFind(_matrix.Ja, dofs[j], a);
                                        _matrix.Elems[a] += local[i, j];
                                    }
                                }
                                _b[dofs[i]] += localB[i];
                            }
                        }
                    }
                }
                // сохранение начала 
                srcxi0 = srcxi;
            }
            srcyi0 = srcyi;
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
