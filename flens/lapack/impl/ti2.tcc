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
       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_TI2_TCC
#define FLENS_LAPACK_IMPL_TI2_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
namespace generic {

template <typename MA>
typename TrMatrix<MA>::IndexType
ti2_impl(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::ElementType ElementType;
    typedef typename TrMatrix<MA>::IndexType   IndexType;

    const Underscore<IndexType>  _;

    const bool upper  = (A.upLo()==Upper);
    const bool noUnit = (A.diag()==NonUnit);
    const IndexType n = A.dim();

    const ElementType  One(1);

    IndexType info = 0;

    if (upper) {
//
//      Compute inverse of upper triangular matrix.
//
        for (IndexType j=1; j<=n; ++j) {
            ElementType Ajj;
            if (noUnit) {
                A(j,j) = One / A(j,j);
                Ajj = -A(j,j);
            } else {
                Ajj = -One;
            }
//
//          Compute elements 1:j-1 of j-th column.
//
            const auto range = _(1,j-1);

            const auto U11 = (noUnit) ? A(range,range).upper()
                                      : A(range,range).upperUnit();
            auto u12 = A(range,j);

            blas::mv(NoTrans, U11, u12);
            u12 *= Ajj;
        }
    } else {
//
//      Compute inverse of lower triangular matrix.
//
        for (IndexType j=n; j>=1; --j) {
            ElementType Ajj;
            if (noUnit) {
                A(j,j) = One / A(j,j);
                Ajj = -A(j,j);
            } else {
                Ajj = -One;
            }
            if (j<n) {
//
//              Compute elements j+1:n of j-th column.
//
                const auto range = _(j+1,n);

                const auto L22 = (noUnit) ? A(range,range).lower()
                                          : A(range,range).lowerUnit();
                auto l21 = A(range,j);

                blas::mv(NoTrans, L22, l21);
                l21 *= Ajj;
            }
        }
    }

    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA>
typename TrMatrix<MA>::IndexType
ti2_impl(TrMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::trti2<IndexType>(getF77Char(A.upLo()),
                                                 getF77Char(A.diag()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
typename TrMatrix<MA>::IndexType
ti2(TrMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename TrMatrix<MA>::NoView   A_org = A;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::ti2_impl(A);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename TrMatrix<MA>::NoView   A_generic = A;

    A = A_org;

    const IndexType _info = external::ti2_impl(A);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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

//-- forwarding ----------------------------------------------------------------

template <typename MA>
typename MA::IndexType
ti2(MA &&A)
{
    typedef typename MA::IndexType  IndexType;

    CHECKPOINT_ENTER;
    IndexType info = ti2(A);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TI2_TCC
