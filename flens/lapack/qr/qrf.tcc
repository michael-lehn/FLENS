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
      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_QR_QRF_TCC
#define FLENS_LAPACK_QR_QRF_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VTAU, typename VWORK>
void
qrf_generic(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    using std::max;
    using std::min;

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = min(m, n);

    ASSERT(tau.length()>=k);
    ASSERT((work.length()>=n) || (work.length()==0));

//
//  Perform and apply workspace query
//
    IndexType nb = ilaenv<T>(1, "GEQRF", "", m, n);

    const IndexType lWorkOpt = n*nb;
    if (work.length()==0) {
        work.resize(max(lWorkOpt,IndexType(1)));
        work(1)=lWorkOpt;
    }

//
//  Quick return if possible
//
    if (k==0) {
        work(1) = 1;
        return;
    }

    IndexType nbMin = 2;
    IndexType nx = 0;
    IndexType iws = n;
    IndexType ldWork = -1;

    if ((nb>1) && (nb<k)) {
//
//      Determine when to cross over from blocked to unblocked code.
//
        nx = max(IndexType(0), IndexType(ilaenv<T>(3, "GEQRF", "", m, n)));
        if (nx<k) {
//
//          Determine if workspace is large enough for blocked code.
//
            ldWork = n;
            iws = ldWork*nb;
            if (work.length()<iws) {
//
//              Not enough workspace to use optimal NB:  reduce NB and
//              determine the minimum value of NB.
//
                nb = work.length() / ldWork;
                nbMin = max(IndexType(2),
                            IndexType(ilaenv<T>(2, "GEQRF", 0, m, n)));
            }
        }
    }

    IndexType i;
    if ((nb>=nbMin) && (nb<k) && (nx<k)) {
        typename GeMatrix<MA>::View Work(n, nb, work);
//
//      Use blocked code initially
//
        for (i=1; i<=k-nx; i+=nb) {
            const IndexType ib = min(k-i+1, nb);
//
//          Compute the QR factorization of the current block
//          A(i:m,i:i+ib-1)
//
            qr2(A(_(i,m),_(i,i+ib-1)), tau(_(i,i+ib-1)), work);
            if (i+ib<=n) {
//
//              Form the triangular factor of the block reflector
//              H = H(i) H(i+1) . . . H(i+ib-1)
//
                auto Tr = Work(_(1,ib),_(1,ib)).upper();
                larft(Forward, ColumnWise,
                      m-i+1,
                      A(_(i,m),_(i,i+ib-1)),
                      tau(_(i,i+ib-1)),
                      Tr);
//
//              Apply H' to A(i:m,i+ib:n) from the left
//
                larfb(Left, Trans, Forward, ColumnWise,
                      A(_(i,m),_(i,i+ib-1)),
                      Tr,
                      A(_(i,m),_(i+ib,n)),
                      Work(_(ib+1,n),_(1,ib)));
            }
        }
    } else {
        i = 1;
    }

//
//  Use unblocked code to factor the last or only block.
//
    if (i<=k) {
        qr2(A(_(i,m),_(i,n)), tau(_(i,k)), work);
    }
    work(1) = iws;
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VTAU, typename VWORK>
void
qrf_native(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const INTEGER    M = A.numRows();
    const INTEGER    N = A.numCols();
    const INTEGER    LDA = A.leadingDimension();
    INTEGER          LWORK = work.length();
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dgeqrf)(&M, &N, A.data(), &LDA,
                            tau.data(),
                            work.data(), &LWORK,
                            &INFO);
        assert(INFO==0);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VTAU, typename VWORK>
void
qrf(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       _A      = A;
    typename DenseVector<VTAU>::NoView  _tau    = tau;
    typename DenseVector<VTAU>::NoView  _work   = work;
#   endif

//
//  Call implementation
//
    qrf_generic(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (work.length()!=_work.length()) {
        // work got resized be the generic implementation
        _work.resize(work.length());
    }
    qrf_native(_A, _tau, _work);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau, _tau, " tau", "_tau")) {
        std::cerr << "CXXLAPACK:  tau = " << tau << std::endl;
        std::cerr << "F77LAPACK: _tau = " << _tau << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename VWORK>
void
qrf(MA &&A, VTAU &&tau, VWORK &&work)
{
    return qrf(A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_QRF_TCC
