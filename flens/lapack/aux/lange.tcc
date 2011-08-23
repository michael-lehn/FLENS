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

/* Based on
 *
       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_AUX_LANGE_TCC
#define FLENS_LAPACK_AUX_LANGE_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

using namespace LANGE;
using std::abs;
using std::max;
using std::min;

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VWORK>
typename MA::ElementType
lange(LANGE::Norm norm, const MA &A, VWORK &&work)
{
    return lange(norm, A, work);
}

//-- lange ---------------------------------------------------------------------
template <typename MA, typename VWORK>
typename GeMatrix<MA>::ElementType
lange(LANGE::Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work)
{
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((norm!=InfinityNorm) || (work>=A.numRows()));

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    if (min(m,n)==T(0)) {
        return 0;
    } else if (norm==MaximumNorm) {
//
//      Find max(abs(A(i,j))).
//
        T value = 0;
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                value = max(value, abs(A(i,j)));
            }
        }
    } else if (norm==OneNorm) {
//
//      Find norm1(A).
//
        T value = 0;

        for (IndexType j=1; j<=n; ++j) {
            T sum = 0;
            for (IndexType i=1; i<=m; ++i) {
                sum += abs(A(i,j));
            }
            value = max(value, sum);
        }
    } else if (norm==InfinityNorm) {
//
//      Find normI(A).
//
        for (IndexType i=1; i<=m; ++i) {
            work(i) = 0;
        }
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=m; ++i) {
                work(i) += abs(A(i,j));
            }
        }
        T value = 0;
        for (IndexType i=1; i<=m; ++i) {
            value = max(value, work(i));
        }
    } else if (norm==FrobeniusNorm) {
//
//      Find normF(A).
//
        T scale = 0;
        T sum = 1;
        for (IndexType j=1; j<=n; ++j) {
            lassq(A(_,j), scale, sum);
        }
    }
    ASSERT(0);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LANGE_TCC
