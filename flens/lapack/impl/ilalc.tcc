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
      INTEGER FUNCTION ILA?LC( M, N, A, LDA )
 *
 *  -- LAPACK auxiliary routine (version 3.2.2)                        --
 *
 *  -- June 2010                                                       --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 */

#ifndef FLENS_LAPACK_IMPL_ILALC_TCC
#define FLENS_LAPACK_IMPL_ILALC_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalc_impl(const GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//  Quick test for the common case where one corner is non-zero.
//
    if (n==0) {
        return n;
    } else if ((A(1,n)!=T(0)) || A(m,n)!=T(0)) {
        return n;
    } else {
//
//      Now scan each column from the end, returning with the first non-zero.
//
        for (IndexType j=n; j>=1; --j) {
            for (IndexType i=1; i<=m; ++i) {
                if (A(i,j)!=T(0)) {
                    return j;
                }
            }
        }
    }
    return 0;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalc_impl(const GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::ilalc<IndexType>(A.numRows(), A.numCols(),
                                       A.data(),
                                       A.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalc(const GeMatrix<MA> &A)
{
    LAPACK_DEBUG_OUT("ilalc");

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    return LAPACK_SELECT::ilalc_impl(A);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ILALC_TCC
