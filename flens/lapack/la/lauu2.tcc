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
       SUBROUTINE DLAUU2( UPLO, N, A, LDA, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011
 */

#ifndef FLENS_LAPACK_LA_LAUU2_TCC
#define FLENS_LAPACK_LA_LAUU2_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- lauu2 [real variant] ------------------------------------------------------
template <typename MA>
typename RestrictTo<IsReal<typename TrMatrix<MA>::ElementType>::value,
                    void>::Type
lauu2_impl(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::ElementType   ElementType;
    typedef typename TrMatrix<MA>::IndexType     IndexType;

    const ElementType            One(1);
    const IndexType              n = A.dim();
    const Underscore<IndexType>  _;
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    if (A.upLo()==Upper) {
//
//      Compute the product U * U**T.
//
        for (IndexType i=1; i<=n; ++i) {
            const ElementType a22 = A(i,i);
            if (i<n) {
                A(i,i) = A(i,_(i,n)) * A(i,_(i,n));

                const auto range1 = _(1,i-1);
                const auto range2 = i;
                const auto range3 = _(i+1,n);

                auto a12       = A(range1,range2);
                const auto A13 = A(range1,range3);
                const auto a23 = A(range2,range3);
                blas::mv(NoTrans, One, A13, a23, a22, a12);
            } else {
                A(_,n) *= a22;
            }
        }
    } else {
//
//      Compute the product L**T * L.
//
        for (IndexType i=1; i<=n; ++i) {
            const ElementType a22 = A(i,i);
            if (i<n) {
                A(i,i) = A(_(i,n),i) * A(_(i,n),i);

                const auto range1 = _(1,i-1);
                const auto range2 = i;
                const auto range3 = _(i+1,n);

                auto a21       = A(range2,range1);
                const auto A31 = A(range3,range1);
                const auto a32 = A(range3,range2);

                blas::mv(Trans, One, A31, a32, a22, a21);
            } else {
                A(n,_) *= a22;
            }
        }
    }
}

//-- lauu2 [complex variant] ------------------------------------------------------
template <typename MA>
typename RestrictTo<IsComplex<typename TrMatrix<MA>::ElementType>::value,
                    void>::Type
lauu2_impl(TrMatrix<MA> &A)
{
    using std::real;
    
    typedef typename TrMatrix<MA>::ElementType   ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType PrimitiveType;
    typedef typename TrMatrix<MA>::IndexType     IndexType;

    const ElementType            COne(1);
    const PrimitiveType          One(1);
    const IndexType              n = A.dim();
    const Underscore<IndexType>  _;
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    if (A.upLo()==Upper) {
//
//      Compute the product U * U**H.
//
        for (IndexType i=1; i<=n; ++i) {
            const PrimitiveType a22 = real(A(i,i));
            if (i<n) {
                A(i,i) = a22*a22 + real(blas::dotc(A(i,_(i+1,n)), A(i,_(i+1,n))));
                imag(A(i,_(i+1,n))) *= -One;
                
                const auto range1 = _(1,i-1);
                const auto range2 = i;
                const auto range3 = _(i+1,n);

                auto a12       = A(range1,range2);
                const auto A13 = A(range1,range3);
                const auto a23 = A(range2,range3);
                blas::mv(NoTrans, COne, A13, a23, ElementType(a22,0), a12);
                imag(A(i,_(i+1,n))) *= -One;
            } else {
                A(_,n) *= a22;
            }
        }
    } else {
//
//      Compute the product L**H * L.
//
        for (IndexType i=1; i<=n; ++i) {
            const PrimitiveType a22 = real(A(i,i));
            if (i<n) {
                A(i,i) = a22*a22 + real(blas::dotc(A(_(i+1,n),i), A(_(i+1,n),i)));
                imag(A(i,_(1,i-1))) *= -One;
                
                const auto range1 = _(1,i-1);
                const auto range2 = i;
                const auto range3 = _(i+1,n);

                auto a21       = A(range2,range1);
                const auto A31 = A(range3,range1);
                const auto a32 = A(range3,range2);

                blas::mv(ConjTrans, COne, A31, a32, ElementType(a22,0), a21);
                imag(A(i,_(1,i-1))) *= -One;
            } else {
                A(n,_) *= a22;
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
lauu2_impl(TrMatrix<MA> &A)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::lauu2<IndexType>(getF77Char(A.upLo()),
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
lauu2(TrMatrix<MA> &A)
{

//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename TrMatrix<MA>::NoView       A_org      = A;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lauu2_impl(A);

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
    external::lauu2_impl(A);

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
lauu2(MA &&A)
{
    lauu2(A);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAUU2_TCC
