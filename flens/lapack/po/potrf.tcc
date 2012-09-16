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
       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
       SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_PO_POTRF_TCC
#define FLENS_LAPACK_PO_POTRF_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- potrf [real variant] ------------------------------------------------------

template <typename MA>
typename SyMatrix<MA>::IndexType
potrf_impl(SyMatrix<MA> &A)
{
    typedef typename SyMatrix<MA>::ElementType  T;
    typedef typename SyMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType n = A.dim();
    const bool upper = (A.upLo()==Upper);

    const T One(1);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Determine the block size for this environment.
//
    const char *upLo = (upper) ? "U" : "L";
    IndexType nb = ilaenv<T>(1, "POGTRF", upLo, n);

    if ((nb<=1) || (nb>=n)) {
//
//      Use unblocked code.
//
        info = potf2(A);
    } else {
//
//      Use blocked code.
//
        if (upper) {
//
//          Compute the Cholesky factorization A = U**T*U.
//
            for (IndexType j=1; j<=n; j+=nb) {
//
//              Partition matrix
//
                IndexType jb = min(nb, n-j+1);

                const auto range1 = _(1,j-1);
                const auto range2 = _(j,j+jb-1);
                const auto range3 = _(j+jb,n);

                const auto A12 = A(range1,range2);
                const auto A13 = A(range1,range3);
                const auto U22 = A(range2,range2).upper();
                auto A22 = A(range2,range2).upper().symmetric();
                auto A23 = A(range2,range3);
//
//              Update and factorize the current diagonal block and test
//              for non-positive-definiteness.
//
                blas::rk(Trans, -One, A12, One, A22);
                info = potf2(A22);
                if (info!=0) {
                    info += j-1;
                    break;
                }
                if (j+jb<=n) {
//
//                  Compute the current block row.
//
                    blas::mm(Trans, NoTrans, -One, A12, A13, One, A23);
                    blas::sm(Left, Trans, One, U22, A23);
                }
            }
        } else {
//
//          Compute the Cholesky factorization A = L*L**T.
//
            for (IndexType j=1; j<=n; j+=nb) {
//
//              Partition matrix
//
                IndexType jb = min(nb, n-j+1);

                const auto range1 = _(1,j-1);
                const auto range2 = _(j,j+jb-1);
                const auto range3 = _(j+jb,n);

                const auto A21 = A(range2,range1);
                const auto A31 = A(range3,range1);
                const auto L22 = A(range2,range2).lower();
                auto A22 = A(range2,range2).lower().symmetric();
                auto A32 = A(range3,range2);
//
//              Update and factorize the current diagonal block and test
//              for non-positive-definiteness.
//
                blas::rk(NoTrans, -One, A21, One, A22);
                info = potf2(A22);
                if (info!=0) {
                    info += j-1;
                    break;
                }
                if (j+jb<=n) {
//
//                  Compute the current block column.
//
                    blas::mm(NoTrans, Trans, -One, A31, A21, One, A32);
                    blas::sm(Right, Trans, One, L22, A32);
                }
            }
        }
    }
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- potrf [real variant] ------------------------------------------------------

template <typename MA>
typename SyMatrix<MA>::IndexType
potrf_impl(SyMatrix<MA> &A)
{
    typedef typename SyMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::potrf<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- potrf [complex variant] ---------------------------------------------------

template <typename MA>
typename HeMatrix<MA>::IndexType
potrf_impl(HeMatrix<MA> &A)
{
    typedef typename HeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::potrf<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}


} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- potrf [real variant] ------------------------------------------------------

template <typename MA>
typename RestrictTo<IsRealSyMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
potrf(MA &&A)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView       _A      = A;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::potrf_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const IndexType _info = external::potrf_impl(_A);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }

#   endif

    return info;
}

#ifdef USE_CXXLAPACK

//-- potrf [complex variant] ---------------------------------------------------

template <typename MA>
typename RestrictTo<IsHeMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
potrf(MA &&A)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

//
//  Call implementation
//
    const IndexType info = external::potrf_impl(A);

    return info;
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_PO_POTRF_TCC
