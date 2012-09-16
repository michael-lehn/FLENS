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
       SUBROUTINE DLAUUM( UPLO, N, A, LDA, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011
 */

#ifndef FLENS_LAPACK_LA_LAUUM_TCC
#define FLENS_LAPACK_LA_LAUUM_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA>
void
lauum_impl(TrMatrix<MA> &A)
{
    using std::min;

    typedef typename TrMatrix<MA>::ElementType   ElementType;
    typedef typename TrMatrix<MA>::IndexType     IndexType;

    const ElementType            One(1);
    const IndexType              n = A.dim();
    const Underscore<IndexType>  _;

    const bool upper = (A.upLo()==Upper);
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
//
//  Determine the block size for this environment.
//
    const char *upLo = (upper) ? "U" : "L";
    const IndexType nb = ilaenv<ElementType>(1, "LAUUM", upLo, n);

    if (nb<=1 || nb>=n) {
//
//      Use unblocked code
//
        lauu2(A);
    } else {
//
//      Use blocked code
//
        if (upper) {
//
//          Compute the product U * U**T.
//
            for (IndexType i=1; i<=n; i+=nb) {
                const IndexType ib = min(nb, n-i+1);

                const auto range1 = _(1,i-1);
                const auto range2 = _(i,i+ib-1);
                const auto range3 = _(i+ib,n);

                auto A12       = A(range1,range2);
                const auto A13 = A(range1,range3);
                auto U22       = A(range2,range2).upper();
                auto A22       = U22.symmetric();
                const auto A23 = A(range2,range3);

                blas::mm(Right, Trans, One, U22, A12);
                lauu2(U22);
                if (i+ib<=n) {
                    blas::mm(NoTrans, Trans, One, A13, A23, One, A12);
                    blas::rk(NoTrans, One, A23, One, A22);
                }
            }
        } else {
//
//          Compute the product L**T * L.
//
            for (IndexType i=1; i<=n; i+=nb) {
                const IndexType ib = min(nb, n-i+1);

                const auto range1 = _(1,i-1);
                const auto range2 = _(i,i+ib-1);
                const auto range3 = _(i+ib,n);

                auto A21       = A(range2,range1);
                const auto A31 = A(range3,range1);
                auto L22       = A(range2,range2).lower();
                auto A22       = L22.symmetric();
                const auto A32 = A(range3,range2);

                blas::mm(Left, Trans, One, L22, A21);
                lauu2(L22);
                if (i+ib<=n) {
                    blas::mm(Trans, NoTrans, One, A32, A31, One, A21);
                    blas::rk(Trans, One, A32, One, A22);
                }
            }
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA>
void
lauum_impl(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::lauum<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
void
lauum(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::IndexType    IndexType;

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.diag()==NonUnit);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename TrMatrix<MA>::NoView       A_org      = A;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lauum_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of generic results
//
    typename TrMatrix<MA>::NoView       A_generic      = A;
//
//  Restore output arguments
//
    A = A_org;

//
//  Compare results
//
    external::lauum_impl(A);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "_A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }

#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA>
void
lauum(MA &&A)
{
    lauum(A);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAUUM_TCC
