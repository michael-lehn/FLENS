/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
       SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_TZRZF_TCC
#define FLENS_LAPACK_IMPL_TZRZF_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- tzrzf [real variant] ------------------------------------------------------

template <typename MA, typename VTAU, typename VWORK>
void
tzrzf_impl(GeMatrix<MA>              &A,
           DenseVector<VTAU>         &tau,
           DenseVector<VWORK>        &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    typedef GeMatrixView<ElementType>           GeView;

    const Underscore<IndexType>  _;

    const ElementType  Zero(0);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    IndexType  lWorkMin, lWorkOpt, nb;

    if (m==0 || m==n) {
        lWorkOpt = 1;
        lWorkMin = 1;
    } else {
//
//      Determine the block size.
//
        nb = ilaenv<ElementType>(1, "GERQF", "", m, n);
        lWorkOpt = m*nb;
        lWorkMin = max(IndexType(1), m);
    }
    if (work.length()==0) {
        work.resize(lWorkOpt);
    }
    work(1) = lWorkOpt;
    ASSERT(work.length()>=lWorkMin);

//
//  Quick return if possible
//
    if (m==0) {
        return;
    } else if (m==n) {
        tau = Zero;
        return;
    }

    IndexType nbMin  = 2;
    IndexType nx     = 1;
    IndexType iws    = m;
    IndexType ldWork = -1;

    if (nb>1 && nb<m) {
//
//      Determine when to cross over from blocked to unblocked code.
//
        nx = max(IndexType(0), ilaenv<ElementType>(3, "GERQF", "", m, n));
        if (nx<m) {
//
//          Determine if workspace is large enough for blocked code.
//
            ldWork = m;
            iws = ldWork*nb;
            if (work.length()<iws) {
//
//              Not enough workspace to use optimal NB:  reduce NB and
//              determine the minimum value of NB.
//
                nb = work.length()/ldWork;
                nbMin = max(IndexType(2),
                            ilaenv<ElementType>(2, "GERQF", "", m, n));
            }
        }
    }

    IndexType mu = 0;
    if (nb>=nbMin && nb<m && nx<m) {
//
//      Use blocked code initially.
//      The last kk rows are handled by the block method.
//
        IndexType m1 = min(m+1, n);
        IndexType ki = ((m-nx-1)/nb)*nb;
        IndexType kk = min(m, ki+nb);

        IndexType i;
        for (i=m-kk+ki+1; i>=m-kk+1; i-=nb) {
            const IndexType ib = min(m-i+1, nb);

            GeView Work(ldWork, ib, work, ldWork);
            auto Tr = Work(_(1,ib),_(1,ib)).lower();
//
//          Compute the TZ factorization of the current block
//          A(i:i+ib-1,i:n)
//
            auto _tau  = tau(_(i,i+ib-1));
            auto _work = work(_(1,ib));
            latrz(n-m, A(_(i,i+ib-1),_(i,n)), _tau, _work);

            if (i>1) {
//
//              Form the triangular factor of the block reflector
//              H = H(i+ib-1) . . . H(i+1) H(i)
//
                larzt(Backward, RowWise, A(_(i,i+ib-1),_(m1,n)), _tau, Tr);
//
//              Apply H to A(1:i-1,i:n) from the right
//
                larzb(Right, NoTrans, Backward, RowWise,
                      A(_(i,i+ib-1),_(m1,n)), Tr, A(_(1,i-1),_(i,n)),
                      Work(_(ib+1,ib+i-1),_));
            }
        }
        mu = i + nb - 1;
    } else {
        mu = m;
    }
//
//  Use unblocked code to factor the last or only block
//
    if (mu>0) {
        latrz(n-m, A(_(1,mu),_), tau(_(1,mu)), work(_(1,mu)));
    }

    work(1) = lWorkOpt;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- tzrzf [real  and complex variant] -----------------------------------------

template <typename MA, typename VTAU, typename VWORK>
void
tzrzf_impl(GeMatrix<MA>              &A,
           DenseVector<VTAU>         &tau,
           DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::tzrzf<IndexType>(A.numRows(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- tzrzf [real and complex variant] ------------------------------------------

template <typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsDenseVector<VTAU>::value
                 && IsDenseVector<VWORK>::value,
         void>::Type
tzrzf(MA           &&A,
      VTAU         &&tau,
      VWORK        &&work)
{
    using std::max;
    using std::min;

    LAPACK_DEBUG_OUT("tzrzf [real]");

//
//  Test the input parameters
//

//
//  Remove references from rvalue types
//
#   if defined(CHECK_CXXLAPACK) || !defined(NDEBUG)

    typedef typename RemoveRef<MA>::Type    MatrixA;
    
#   endif

#   ifndef NDEBUG

    typedef typename MatrixA::IndexType     IndexType;
    
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(n>=m);
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(tau.length()==m);

    if (work.length()>0) {
        ASSERT(work.length()>=max(IndexType(1),m));
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<VTAU>::Type  VectorTau;
    typedef typename RemoveRef<VWORK>::Type VectorWork;
    
    typename MatrixA::NoView        A_org    = A;
    typename VectorTau::NoView      tau_org  = tau;
    typename VectorWork::NoView     work_org = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::tzrzf_impl(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView        A_generic    = A;
    typename VectorTau::NoView      tau_generic  = tau;
    typename VectorWork::NoView     work_generic = work;

//
//  restore output arguments
//
    A    = A_org;
    tau  = tau_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::tzrzf_impl(A, tau, work);

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
        std::cerr << "error in: tzrzf.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: tzrzf.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TZRZF_TCC
