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

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = min(m, n);
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
                            IndexType(ilaenv<T>(2, "GEQRF", "", m, n)));
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
            qr2(A(_(i,m),_(i,i+ib-1)), tau(_(i,i+ib-1)), work(_(1,ib)));
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
        qr2(A(_(i,m),_(i,n)), tau(_(i,k)), work(_(1,n-i+1)));
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
    using std::min;
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = min(m,n);

    ASSERT(tau.length()==k);
    ASSERT(work.length()>=n || work.length()==IndexType(0));
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename DenseVector<VTAU>::NoView  tau_org    = tau;
    typename DenseVector<VWORK>::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    qrf_generic(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename GeMatrix<MA>::NoView       A_generic      = A;
    typename DenseVector<VTAU>::NoView  tau_generic    = tau;
    typename DenseVector<VWORK>::NoView work_generic   = work;

    A    = A_org;
    tau  = tau_org;

    // if the generic implementation resized work due to a work size query
    // we must not restore the work array
    if (work_org.length()>0) {
        work = work_org;
    } else {
        work = 0;
    }
//
//  Compare results
//
    qrf_native(A, tau, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename VWORK>
void
qrf(MA &&A, VTAU &&tau, VWORK &&work)
{
    CHECKPOINT_ENTER;
    qrf(A, tau, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_QRF_TCC
