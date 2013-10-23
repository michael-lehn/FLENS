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
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Based on
 *
      INTEGER FUNCTION ILA?LR( M, N, A, LDA )
 *
 *  -- LAPACK auxiliary routine (version 3.2.2)                        --
 *
 *  -- June 2010                                                       --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 */

#ifndef FLENS_LAPACK_LA_ILALR_TCC
#define FLENS_LAPACK_LA_ILALR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalr_impl(const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::IndexType    IndexType;
    typedef typename GeMatrix<MA>::ElementType  T;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//  Quick test for the common case where one corner is non-zero.
//
    if (m==0) {
        return m;
    }
    if ((A(m,1)!=T(0)) || A(m,n)!=T(0)) {
        return m;
    }
//
//  Scan up each column tracking the last zero row seen.
//
    IndexType lastRow = 0;
    for (IndexType j=1; j<=n; ++j) {
        for (IndexType i=m; i>=1; --i) {
            if (A(i,j)!=T(0)) {
                lastRow = max(lastRow, i);
                break;
            }
        }
    }
    return lastRow;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalr_impl(const GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::ilalr<IndexType>(A.numRows(), A.numCols(),
                                       A.data(),
                                       A.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
typename GeMatrix<MA>::IndexType
ilalr(const GeMatrix<MA> &A)
{
    LAPACK_DEBUG_OUT("ilalr");

    typedef typename GeMatrix<MA>::IndexType   IndexType;

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType info = LAPACK_SELECT::ilalr_impl(A);

#   ifdef CHECK_CXXLAPACK
    const IndexType _info = external::ilalr_impl(A);

    if (info!=_info) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;

        std::cerr << "A = " << A << std::endl;

        ASSERT(0);
    }
#   endif

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_ILALR_TCC
