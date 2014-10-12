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

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

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

                const auto A_  = A(_(1,i),_(i+1,n));
                auto       w_  = W(i,_(iw+1,nb));
                auto       A_i = A(_(1,i),i);

                blas::conj(w_);
                blas::mv(NoTrans, -One, A_, w_, One, A_i);
                blas::conj(w_);

                const auto W_  = W(_(1,i),_(iw+1,nb));
                auto       Ai_ = A(i,_(i+1,n));

                blas::conj(Ai_);
                blas::mv(NoTrans, -One, W_, Ai_, One, A_i);
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
                const auto A11_  = A(_(1,i-1),_(1,i-1)).upper().hermitian();
                const auto Ai_   = A(_(1,i-1),i);
                auto       W1iw_ = W(_(1,i-1),iw);

                blas::mv(One, A11_, Ai_, Zero, W1iw_);
                if (i<n) {
                    const auto W_    = W(_(1,i-1),_(iw+1,nb));
                    auto       W2iw_ = W(_(i+1,n),iw);

                    blas::mv(ConjTrans, One, W_, Ai_, Zero, W2iw_);

                    const auto A12_ = A(_(1,i-1),_(i+1,n));

                    blas::mv(NoTrans, -One, A12_, W2iw_, One, W1iw_);
                    blas::mv(ConjTrans, One, A12_, Ai_, Zero, W2iw_);
                    blas::mv(NoTrans, -One, W_, W2iw_, One, W1iw_);
                }
                W1iw_ *= tau(i-1);
                alpha = -Half * tau(i-1) * blas::dotc(W1iw_, Ai_);
                W1iw_ += alpha*Ai_;
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

            const auto A_  = A(_(i,n),_(1,i-1));
            auto       w_  = W(i,_(1,i-1));
            auto       A_i = A(_(i,n),i);

            blas::conj(w_);
            blas::mv(NoTrans, -One, A_, w_, One, A_i);
            blas::conj(w_);

            const auto W_  = W(_(i,n),_(1,i-1));
            auto       Ai_ = A(i,_(1,i-1));

            blas::conj(Ai_);
            blas::mv(NoTrans, -One, W_, Ai_, One, A_i);
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
                const auto A22_ = A(_(i+1,n),_(i+1,n)).lower().hermitian();
                const auto Ai_  = A(_(i+1,n),i);
                auto       W2i_ = W(_(i+1,n),i);

                blas::mv(One, A22_, Ai_, Zero, W2i_);

                const auto W_   = W(_(i+1,n),_(1,i-1));
                auto       W1i_ = W(_(1,i-1),i);

                blas::mv(ConjTrans, One, W_, Ai_, Zero, W1i_);

                const auto A21_ = A(_(i+1,n),_(1,i-1));

                blas::mv(NoTrans, -One, A21_, W1i_, One, W2i_);
                blas::mv(ConjTrans, One, A21_, Ai_, Zero, W1i_);
                blas::mv(NoTrans, -One, W_, W1i_, One, W2i_);
                W2i_ *= tau(i);
                alpha = -Half * tau(i) * blas::dotc(W2i_, Ai_);
                W2i_ += alpha*Ai_;
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

