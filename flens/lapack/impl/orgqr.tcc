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
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_ORGQR_TCC
#define FLENS_LAPACK_IMPL_ORGQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orgqr_impl(IndexType                 k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;
    const T  Zero(0);

    const Underscore<IndexType> _;
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
//
//  Perform and apply workspace query
//
    IndexType nb = ilaenv<T>(1, "ORGQR", "", m, n, k);
    const IndexType lWorkOpt = max(IndexType(1), n) * nb;
    if (work.length()==0) {
        work.resize(max(lWorkOpt,IndexType(1)));
        work(1)=lWorkOpt;
    }
//
//  Quick return if possible
//
    if (n<=0) {
        work(1) = 1;
        return;
    }

    IndexType nbMin = 2;
    IndexType nx = 0;
    IndexType iws = n;
    IndexType ldWork = 1;

    if ((nb>1) && (nb<k)) {
//
//      Determine when to cross over from blocked to unblocked code.
//
        nx = max(IndexType(0), IndexType(ilaenv<T>(3, "ORGQR", "", m, n, k)));
        if (nx<k) {
//
//          Determine if workspace is large enough for blocked code.
//
            ldWork = n;
            iws = ldWork *nb;
            if (work.length()<iws) {
//
//              Not enough workspace to use optimal NB:  reduce NB and
//              determine the minimum value of NB.
//
                nb = work.length() / ldWork;
                nbMin = max(IndexType(2),
                            IndexType(ilaenv<T>(2, "ORGQR", "", m, n, k)));
            }
        }
    }

    IndexType ki = -1,
              kk = -1;

    if ((nb>=nbMin) && (nb<k) && (nx<k)) {
//
//      Use blocked code after the last block.
//      The first kk columns are handled by the block method.
//
        ki = ((k-nx-1)/nb)*nb;
        kk = min(k, ki+nb);
//
//      Set A(1:kk,kk+1:n) to zero.
//
        A(_(1,kk),_(kk+1,n)) = Zero;
    } else {
        kk = 0;
    }

//
//  Use unblocked code for the last or only block.
//
    if (kk<n) {
        org2r(k-kk, A(_(kk+1,m),_(kk+1,n)), tau(_(kk+1, k)), work(_(1,n-kk)));
    }

    if (kk>0) {
        typename GeMatrix<MA>::View Work(n, nb, work);
//
//      Use blocked code
//
        for (IndexType i=ki+1; i>=1; i-=nb) {
            const IndexType ib = min(nb, k-i+1);
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
//              Apply H to A(i:m,i+ib:n) from the left
//
                larfb(Left, NoTrans, Forward, ColumnWise,
                      A(_(i,m),_(i,i+ib-1)),
                      Tr,
                      A(_(i,m),_(i+ib,n)),
                      Work(_(ib+1,n),_(1,ib)));
            }
//
//          Apply H to rows i:m of current block
//
            org2r(ib, A(_(i,m),_(i,i+ib-1)), tau(_(i,i+ib-1)), work(_(1,ib)));
//
//          Set rows 1:i-1 of current block to zero
//
            A(_(1,i-1),_(i,i+ib-1)) = Zero;
        }
    }
    work(1) = iws;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orgqr_impl(IndexType                  k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    cxxlapack::orgqr<IndexType>(A.numRows(),
                                A.numCols(),
                                k,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orgqr(IndexType                 k,
      GeMatrix<MA>              &A,
      const DenseVector<VTAU>   &tau,
      DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));
    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(tau.length()==k);
    ASSERT((work.length()==0) || (work.length()>=A.numCols()));

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(n<=m);
    ASSERT(k<=n);
    ASSERT(0<=k);
#   endif 

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename DenseVector<VWORK>::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::orgqr_impl(k, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename GeMatrix<MA>::NoView       A_generic      = A;
    typename DenseVector<VWORK>::NoView work_generic   = work;

    A = A_org;

    if (work_org.length()!=0) {
        work = work_org;
    } else {
        work = ElementType(0);
    }
//
//  Compare results
//
    external::orgqr_impl(k, A, tau, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: orgqr.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: orgqr.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orgqr(IndexType     k,
      MA            &&A,
      const VTAU    &tau,
      VWORK         &&work)
{
    orgqr(k, A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ORGQR_TCC
