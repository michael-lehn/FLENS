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
       SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HE_TRD_TCC
#define FLENS_LAPACK_HE_TRD_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA>
typename MA::IndexType
trd_wsq_impl(const HeMatrix<MA> &A)
{
    typedef typename MA::IndexType      IndexType;
    typedef typename MA::ElementType    T;

    const IndexType n        = A.dim();
    const char      upLo[1]  = { getF77BlasChar(A.upLo()) };

//
//  Determine the block size.
//
    const IndexType nb       = ilaenv<T>(1, "HETRD", upLo, n);

    return n*nb;
}

template <typename MA, typename VD, typename VE, typename VTAU, typename VWORK>
void
trd_impl(HeMatrix<MA>        &A,
         DenseVector<VD>     &d,
         DenseVector<VE>     &e,
         DenseVector<VTAU>   &tau,
         DenseVector<VWORK>  &work)
{
    using std::real;

    typedef typename MA::IndexType                   IndexType;
    typedef typename MA::ElementType                 T;
    typedef typename ComplexTrait<T>::PrimitiveType  PT;

    const PT One(1);
    const T  COne(1);

    const Underscore<IndexType> _;

    const char      upLo[1] = { getF77BlasChar(A.upLo()) };
    const IndexType n       = A.numRows();

//
//  Determine the block size.
//
    IndexType nb             = ilaenv<T>(1, "HETRD", upLo, n);
    const IndexType lWorkOpt = trd_wsq(A);

    if (work.length()==0) {
        work.resize(lWorkOpt);
    }

    const IndexType lWork = work.length();

    IndexType  ldWork = -1;

//
//  Quick return if possible
//
    if (n==0) {
        work(1) = 1;
        return;
    }

    IndexType nx = n;
    IndexType iws = 1;

    if (nb>1 && nb<n) {
//
//      Determine when to cross over from blocked to unblocked code
//      (last block is always handled by unblocked code).
//
        nx = max(nb, ilaenv<T>(3, "HETRD", upLo, n));

        if (nx<n) {
//
//          Determine if workspace is large enough for blocked code.
//
            ldWork = n;
            iws    = ldWork*nb;

            if (lWork < iws) {
//
//              Not enough workspace to use optimal NB:  determine the
//              minimum value of NB, and reduce NB or force use of
//              unblocked code by setting NX = N.
//
                nb = max(lWork/ldWork, IndexType(1));

                IndexType nbMin = ilaenv<T>(2, "HETRD", upLo, n);

                if (nb<nbMin) {
                    nx = n;
                }
            }
        } else {
            nx = n;
        }
    } else {
        nb = 1;
    }

    if (A.upLo()==Upper) {
//
//      Reduce the upper triangle of A.
//      Columns 1:kk are handled by the unblocked method.
//
        IndexType kk = n - ((n-nx+nb-1)/nb)*nb;

        for (IndexType i=n-nb+1; i>=kk+1; i-=nb) {
//
//          Reduce columns i:i+nb-1 to tridiagonal form and form the
//          matrix W which is needed to update the unreduced part of
//          the matrix
//
            auto A_   = A(_(1,i+nb-1),_(1,i+nb-1)).upper().hermitian();
            auto e_   = e(_(1,i+nb-2));
            auto tau_ = tau(_(1,i+nb-2));

            ASSERT(ldWork!=-1);
            typename GeMatrix<MA>::View  Work(i+nb-1, nb, work, ldWork);

            latrd(A_, e_, tau_, Work);

//
//          Update the unreduced submatrix A(1:i-1,1:i-1), using an
//          update of the form:  A := A - V*W**H - W*V**H
//
            auto __A = A(_(1,i-1),_(1,i-1)).upper().hermitian();
            auto __V = A(_(1,i-1),_(i,i+nb-1));
            auto __W = Work(_(1,i-1),_);

            blas::r2k(NoTrans, -COne, __V, __W, One, __A);

//
//          Copy superdiagonal elements back into A, and diagonal
//          elements into D
//
            for (IndexType j=i; j<=i+nb-1; ++j) {
                A(j-1,j) = e(j-1);
                d(j)     = real(A(j,j));
            }
        }
//
//      Use unblocked code to reduce the last or only block
//
        auto A_   = A(_(1,kk),_(1,kk)).upper().hermitian();
        auto d_   = d(_(1,kk));
        auto e_   = e(_(1,kk-1));
        auto tau_ = tau(_(1,kk-1));

        td2(A_, d_, e_, tau_);
    } else {
//
//      Reduce the lower triangle of A
//
        IndexType i;
        for (i=1; i<=n-nx; i+=nb) {
//
//          Reduce columns i:i+nb-1 to tridiagonal form and form the
//          matrix W which is needed to update the unreduced part of
//          the matrix
//
            auto A_   = A(_(i,n),_(i,n)).lower().hermitian();
            auto e_   = e(_(i,n-1));
            auto tau_ = tau(_(i,n-1));

            ASSERT(ldWork!=-1);
            typename GeMatrix<MA>::View  Work(n-i+1, nb, work, ldWork);

            latrd(A_, e_, tau_, Work);
//
//          Update the unreduced submatrix A(i+nb:n,i+nb:n), using
//          an update of the form:  A := A - V*W**H - W*V**H
//
            auto __A = A(_(i+nb,n),_(i+nb,n)).lower().hermitian();
            auto __V = A(_(i+nb,n),_(i,i+nb-1));

            typename GeMatrix<MA>::View  __Work(n-i-nb+1, nb,
                                                work(_(nb+1, lWork)),
                                                ldWork);

            blas::r2k(NoTrans, -COne, __V, __Work, One, __A);

//
//          Copy subdiagonal elements back into A, and diagonal
//          elements into D
//
            for (IndexType j=i; j<=i+nb-1; ++j) {
                A(j+1, j) = e(j);
                d(j)      = real(A(j,j));
            }
        }
//
//      Use unblocked code to reduce the last or only block
//
        auto A_   = A(_(i,n),_(i,n)).lower().hermitian();
        auto d_   = d(_(i,n));
        auto e_   = e(_(i,n-1));
        auto tau_ = tau(_(i,n-1));

        td2(A_, d_, e_, tau_);
    }

    work(1) = lWorkOpt;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (he)trf [complex variant] -------------------------------------------------

template <typename MA>
typename MA::IndexType
trd_wsq_impl(const HeMatrix<MA> &A)
{
    using std::real;

    typedef typename HeMatrix<MA>::IndexType         IndexType;
    typedef typename HeMatrix<MA>::ElementType       T;
    typedef typename ComplexTrait<T>::PrimitiveType  PT;

    T           DUMMY;
    PT          RDUMMY;
    T           WORK;
    IndexType   LWORK = -1;

    cxxlapack::hetrd<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                &DUMMY, A.leadingDimension(),
                                &RDUMMY,
                                &RDUMMY,
                                &DUMMY,
                                &WORK,
                                LWORK);
    return IndexType(real(WORK));
}

template <typename MA, typename VD, typename VE, typename VTAU, typename VWORK>
void
trd_impl(HeMatrix<MA>        &A,
         DenseVector<VD>     &d,
         DenseVector<VE>     &e,
         DenseVector<VTAU>   &tau,
         DenseVector<VWORK>  &work)
{
    typedef typename HeMatrix<MA>::IndexType   IndexType;

    if (work.length()==0) {
        work.resize(trd_wsq_impl(A));
    }

    cxxlapack::hetrd<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                A.data(), A.leadingDimension(),
                                d.data(),
                                e.data(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (he)trd -------------------------------------------------------------------

template <typename MA, typename VD, typename VE, typename VTAU, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VD>::value
                 && IsRealDenseVector<VE>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
trd(MA      &&A,
    VD      &&d,
    VE      &&e,
    VTAU    &&tau,
    VWORK   &&work)
{
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VD>::Type    VectorD;
    typedef typename RemoveRef<VE>::Type    VectorE;
    typedef typename RemoveRef<VTAU>::Type  VectorTau;
    typedef typename RemoveRef<VWORK>::Type VectorWork;
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
    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==0 || work.length()>=1);
#   endif

//
//  Resize output arguments if they are empty and needed
//
    if (work.length()==0) {
        work.resize(trd_wsq(A));
    }
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView      A_org     = A;
    typename VectorD::NoView      d_org     = d;
    typename VectorE::NoView      e_org     = e;
    typename VectorTau::NoView    tau_org   = tau;
    typename VectorWork::NoView   work_org  = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::trd_impl(A, d, e, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView      A_generic     = A;
    typename VectorD::NoView      d_generic     = d;
    typename VectorE::NoView      e_generic     = e;
    typename VectorTau::NoView    tau_generic   = tau;
    typename VectorWork::NoView   work_generic  = work;

    A     = A_org;
    d     = d_org;
    e     = e_org;
    tau   = tau_org;
    work  = work_org;

    external::trd_impl(A, d, e, tau, work);

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

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- (he)trd [workspace query] -------------------------------------------------

template <typename MA>
typename RestrictTo<IsHeMatrix<MA>::value,
         typename MA::IndexType>::Type
trd_wsq(const MA &A)
{
#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif

    const auto ws = LAPACK_SELECT::trd_wsq_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const auto optWorkSize = external::trd_wsq_impl(A);
    if (! isIdentical(optWorkSize, ws, "optWorkSize", "ws")) {
        ASSERT(0);
    }
#   endif

    return ws;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_HE_TRD_TCC
