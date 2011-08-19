/*
 *   Copyright (c) 2011, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *   
 *   $COPYRIGHT$
 *   
 *   Additional copyrights may follow
 *   
 *   $HEADER$
 *   
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *   
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *     
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *     
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *     
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */

#ifndef FLENS_LAPACK_GESV_TF2_TCC
#define FLENS_LAPACK_GESV_TF2_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VP>
typename MA::IndexType
tf2(MA &&A, VP &&piv)
{
    return tf2(A, piv);
}

//-- getf2 ---------------------------------------------------------------------
template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
tf2(GeMatrix<MA> &A, DenseVector<VP> &piv)
{
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.inc()==1);

    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    typedef Range<IndexType>    Range;
    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    IndexType info = 0;

//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return info;
    }
//
//     Compute machine safe minimum 
//
    T safeMin = lamch<T>(SafeMin);

    for (IndexType j=1; j<=std::min(m,n); ++j) {
//
//      Row range of current submatrix A(j:M, j:N)
//
        const Range rows(j, m);
//
//      Row and column range of trailing submatrix A(j+1:M, j+1:N)
//
        const Range _rows(j+1, m);
        const Range _cols(j+1, n);
//
//      Find pivot and test for singularity.
//
        IndexType jp = j + blas::iamax(A(rows,j));
        piv(j) = jp;
        if (A(jp, j)!=T(0)) {
//
//          Apply the interchange to columns 1:N.
//
            if (j!=jp) {
                blas::swap(A(j,_), A(jp,_));
            }
//
//          Compute elements J+1:M of J-th column.
//
            if (j<m) {
                if (abs(A(j,j))>=safeMin) {
                    blas::scal(T(1)/A(j, j), A(_rows,j));
                } else {
                    for (IndexType i=1; i<=m-j; ++i) {
                        A(j+i,j) = A(j+i,j)/A(j,j);
                    }
                }
            }
        } else {
            if (info==0) {
                info = j;
            }
        }
//
//      Update trailing submatrix A(j+1:M, j+1:N)
//
        blas::r(T(-1),
                A(_rows,     j),
                A(    j, _cols),
                A(_rows, _cols));
    }

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_TF2_TCC
