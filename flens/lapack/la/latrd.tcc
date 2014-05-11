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
       SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LATRD_TCC
#define FLENS_LAPACK_LA_LATRD_TCC 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VE, typename VTAU, typename MW>
void
latrd_impl(HeMatrix<MA>       &A,
           DenseVector<VE>    &e,
           DenseVector<VTAU>  &tau,
           GeMatrix<MW>       &W)
{
    using std::real;

    typedef typename MA::ElementType  T;
    typedef typename MA::IndexType    IndexType;

    const T Zero(0), Half(0.5), One(1);

    const Underscore<IndexType> _;

    const IndexType n  = A.dim();
    const IndexType nb = W.numCols();
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    if (A.upLo()==Upper) {
//
//      Reduce last NB columns of upper triangle
//
        for (IndexType i=n; i>n-nb; --i) {
            const IndexType iw = i-n+nb;
            if (i<n) {
//
//              Update A(1:i,i)
//
                A(i,i) = real(A(i,i));

                const auto _A  = A(_(1,i),_(i+1,n));
                auto       _w  = W(i,_(iw+1,nb));
                auto       _Ai = A(_(1,i),i);

                blas::conj(_w);
                blas::mv(NoTrans, -One, _A, _w, One, _Ai);
                blas::conj(_w);

                const auto _W  = W(_(1,i),_(iw+1,nb));
                auto       Ai_ = A(i,_(i+1,n));

                blas::conj(Ai_);
                blas::mv(NoTrans, -One, _W, Ai_, One, _Ai);
                blas::conj(Ai_);
                A(i,i) = real(A(i,i));
            }
            if (i>1) {
//
//              Generate elementary reflector H(i) to annihilate
//              A(1:i-2,i)
//
                T alpha = A(i-1,i);
                larfg(i-1, alpha, A(_(1,i-2),i), tau(i-1));
                e(i-1) = real(alpha);
                A(i-1,i) = One;
//
//              Compute W(1:i-1,i)
//
                const auto _A11  = A(_(1,i-1),_(1,i-1)).upper().hermitian();
                const auto _Ai   = A(_(1,i-1),i);
                auto       _W1iw = W(_(1,i-1),iw);

                blas::mv(One, _A11, _Ai, Zero, _W1iw);
                if (i<n) {
                    const auto _W    = W(_(1,i-1),_(iw+1,nb));
                    auto       _W2iw = W(_(i+1,n),iw);

                    blas::mv(ConjTrans, One, _W, _Ai, Zero, _W2iw);

                    const auto _A12 = A(_(1,i-1),_(i+1,n));

                    blas::mv(NoTrans, -One, _A12, _W2iw, One, _W1iw);
                    blas::mv(ConjTrans, One, _A12, _Ai, Zero, _W2iw);
                    blas::mv(NoTrans, -One, _W, _W2iw, One, _W1iw);
                }
                _W1iw *= tau(i-1);
                alpha = -Half * tau(i-1) * blas::dotc(_W1iw, _Ai);
                _W1iw += alpha*_Ai;
            }

        }
    } else {
//
//      Reduce first NB columns of lower triangle
//
        for (IndexType i=1; i<=nb; ++i) {
//
//          Update A(i:n,i)
//
            A(i,i) = real(A(i,i));

            const auto _A  = A(_(i,n),_(1,i-1));
            auto       _w  = W(i,_(1,i-1));
            auto       _Ai = A(_(i,n),i);

            blas::conj(_w);
            blas::mv(NoTrans, -One, _A, _w, One, _Ai);
            blas::conj(_w);

            const auto _W  = W(_(i,n),_(1,i-1));
            auto       Ai_ = A(i,_(1,i-1));

            blas::conj(Ai_);
            blas::mv(NoTrans, -One, _W, Ai_, One, _Ai);
            blas::conj(Ai_);
            A(i,i) = real(A(i,i));

            if (i<n) {
//
//              Generate elementary reflector H(i) to annihilate
//              A(i+2:n,i)
//
                T alpha = A(i+1,i);
                larfg(n-i, alpha, A(_(i+2,n),i), tau(i));
                e(i) = real(alpha);
                A(i+1,i) = One;
//
//              Compute W(i+1:n,i)
//
                const auto _A22 = A(_(i+1,n),_(i+1,n)).lower().hermitian();
                const auto _Ai  = A(_(i+1,n),i);
                auto       _W2i = W(_(i+1,n),i);

                blas::mv(One, _A22, _Ai, Zero, _W2i);

                const auto _W   = W(_(i+1,n),_(1,i-1));
                auto       _W1i = W(_(1,i-1),i);

                blas::mv(ConjTrans, One, _W, _Ai, Zero, _W1i);

                const auto _A21 = A(_(i+1,n),_(1,i-1));

                blas::mv(NoTrans, -One, _A21, _W1i, One, _W2i);
                blas::mv(ConjTrans, One, _A21, _Ai, Zero, _W1i);
                blas::mv(NoTrans, -One, _W, _W1i, One, _W2i);
                _W2i *= tau(i);
                alpha = -Half * tau(i) * blas::dotc(_W2i, _Ai);
                _W2i += alpha*_Ai;
            }
        }
    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VE, typename VTAU, typename MW>
void
latrd_impl(HeMatrix<MA>       &A,
           DenseVector<VE>    &e,
           DenseVector<VTAU>  &tau,
           GeMatrix<MW>       &W)
{
    cxxlapack::latrd(getF77Char(A.upLo()),
                     A.numRows(),
                     W.numCols(),
                     A.data(),
                     A.leadingDimension(),
                     e.data(),
                     tau.data(),
                     W.data(),
                     W.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================
//
//  Complex variant
//
template <typename MA, typename VE, typename VTAU, typename MW>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VE>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexGeMatrix<MW>::value,
         void>::Type
latrd(MA      &&A,
      VE      &&e,
      VTAU    &&tau,
      MW      &&W)
{
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VE>::Type    VectorE;
    typedef typename RemoveRef<VTAU>::Type  VectorTau;
    typedef typename RemoveRef<MW>::Type    MatrixW;
#   endif

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(e.firstIndex()==1);
    ASSERT(e.length()==A.numRows()-1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(tau.length()==A.numRows()-1);
    ASSERT(W.firstRow()==1);
    ASSERT(W.firstCol()==1);
    ASSERT(W.numRows()==A.numRows());
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView      A_org   = A;
    typename VectorE::NoView      e_org   = e;
    typename VectorTau::NoView    tau_org = tau;
    typename MatrixW::NoView      W_org   = W;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::latrd_impl(A, e, tau, W);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView      A_generic   = A;
    typename VectorE::NoView      e_generic   = e;
    typename VectorTau::NoView    tau_generic = tau;
    typename MatrixW::NoView      W_generic   = W;

    A   = A_org;
    e   = e_org;
    tau = tau_org;
    W   = W_org;

    external::latrd_impl(A, e, tau, W);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "char(A.upLo()) = " << char(A.upLo()) << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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

    if (! isIdentical(W_generic, W, "W_generic", "W")) {
        std::cerr << "CXXLAPACK: W_generic = " << W_generic << std::endl;
        std::cerr << "F77LAPACK: W = " << W << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LATRD_TCC

