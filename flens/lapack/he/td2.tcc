/*
 *   Copyright (c) 2014, Michael Lehn
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
       SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HE_TD2_TCC
#define FLENS_LAPACK_HE_TD2_TCC 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VD, typename VE, typename VTAU>
void
td2_impl(HeMatrix<MA>        &A,
         DenseVector<VD>     &d,
         DenseVector<VE>     &e,
         DenseVector<VTAU>   &tau)
{
    using std::real;

    typedef typename MA::ElementType  T;
    typedef typename MA::IndexType    IndexType;

    const Underscore<IndexType> _;

    const T Zero(0), Half(0.5), One(1);

    const IndexType n = A.numRows();

//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    if (A.upLo()==Upper) {
//
//      Reduce the upper triangle of A
//
        A(n,n) = real(A(n,n));
        for (IndexType i=n-1; i>=1; --i) {
//
//          Generate elementary reflector H(i) = I - tau * v * v**H
//          to annihilate A(1:i-1,i+1)
//
            T  alpha = A(i,i+1);

            T taui;
            larfg(i, alpha, A(_(1,i-1),i+1), taui);
            e(i) = real(alpha);

            if (taui!=Zero) {
//
//              Apply H(i) from both sides to A(1:i,1:i)
//
                A(i,i+1) = One;
//
//              Compute  x := tau * A * v  storing x in TAU(1:i)
//
                auto       _A = A(_(1,i),_(1,i)).upper().hermitian();
                const auto _v = A(_(1,i),i+1);
                auto       _x = tau(_(1,i));

                blas::mv(taui, _A, _v, Zero, _x);
//
//              Compute  w := x - 1/2 * tau * (x**H * v) * v
//
                alpha = -Half * taui * blas::dotc(_x, _v);
                _x += alpha*_v;
//
//              Apply the transformation as a rank-2 update:
//                 A := A - v * w**H - w * v**H
//
                blas::r2(-One, _v, _x, _A);

            } else {
                A(i,i) = real(A(i,i));
            }
            A(i,i+1) = e(i);
            d(i+1) = real(A(i+1,i+1));
            tau(i) = taui;
        }
        d(1) = real(A(1,1));
    } else {
//
//      Reduce the lower triangle of A
//
        A(1,1) = real(A(1,1));
        for (IndexType i=1; i<=n-1; ++i) {
//
//          Generate elementary reflector H(i) = I - tau * v * v**H
//          to annihilate A(i+2:n,i)
//
            T alpha =A(i+1,i);

            T taui;
            larfg(n-i, alpha, A(_(i+2,n),i), taui);
            e(i) = real(alpha);

            if (taui!=Zero) {
//
//              Apply H(i) from both sides to A(i+1:n,i+1:n)
//
                A(i+1,i) = One;
//
//              Compute  x := tau * A * v  storing y in TAU(i:n-1)
//
                auto       _A = A(_(i+1,n),_(i+1,n)).lower().hermitian();
                const auto _v = A(_(i+1,n),i);
                auto       _x = tau(_(i,n-1));

                blas::mv(taui, _A, _v, Zero, _x);
//
//              Compute  w := x - 1/2 * tau * (x**H * v) * v
//
                alpha = -Half * taui * blas::dotc(_x, _v);
                _x += alpha*_v;
//
//              Apply the transformation as a rank-2 update:
//                 A := A - v * w**H - w * v**H
//
                blas::r2(-One, _v, _x, _A);

            } else {
                A(i+1,i+1) = real(A(i+1,i+1));
            }
            A(i+1,i) = e(i);
            d(i) = real(A(i,i));
            tau(i) = taui;
        }
        d(n) = real(A(n,n));
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VD, typename VE, typename VTAU>
void
td2_impl(HeMatrix<MA>        &A,
         DenseVector<VD>     &d,
         DenseVector<VE>     &e,
         DenseVector<VTAU>   &tau)
{
    typedef typename HeMatrix<MA>::IndexType   IndexType;

    cxxlapack::hetd2<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                A.data(), A.leadingDimension(),
                                d.data(),
                                e.data(),
                                tau.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VD, typename VE, typename VTAU>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VD>::value
                 && IsRealDenseVector<VE>::value
                 && IsComplexDenseVector<VTAU>::value,
         void>::Type
td2(MA      &&A,
    VD      &&d,
    VE      &&e,
    VTAU    &&tau)
{
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VD>::Type    VectorD;
    typedef typename RemoveRef<VE>::Type    VectorE;
    typedef typename RemoveRef<VTAU>::Type  VectorTau;
#   endif

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(d.firstIndex()==1);
    ASSERT(d.length()==A.numRows());
    ASSERT(e.firstIndex()==1);
    ASSERT(e.length()==A.numRows()-1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(tau.length()==A.numRows()-1);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView      A_org     = A;
    typename VectorD::NoView      d_org     = d;
    typename VectorE::NoView      e_org     = e;
    typename VectorTau::NoView    tau_org   = tau;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::td2_impl(A, d, e, tau);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView      A_generic     = A;
    typename VectorD::NoView      d_generic     = d;
    typename VectorE::NoView      e_generic     = e;
    typename VectorTau::NoView    tau_generic   = tau;

    A     = A_org;
    d     = d_org;
    e     = e_org;
    tau   = tau_org;

    external::td2_impl(A, d, e, tau);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(d_generic, d, "d_generic", "d")) {
        std::cerr << "CXXLAPACK: d_generic = " << d_generic << std::endl;
        std::cerr << "F77LAPACK: d = " << d << std::endl;
        failed = true;
    }

    if (! isIdentical(e_generic, e, "e_generic", "e")) {
        std::cerr << "CXXLAPACK: e_generic = " << e_generic << std::endl;
        std::cerr << "F77LAPACK: e = " << e << std::endl;
        failed = true;
    }

    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "A_org.upLo() = " << char(A_org.upLo()) << std::endl;
        std::cerr << "d_org = " << d_org << std::endl;
        std::cerr << "e_org = " << e_org << std::endl;
        std::cerr << "tau_org = " << tau_org << std::endl;
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HE_TD2_TCC
