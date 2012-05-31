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
       SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_ORGLQ_TCC
#define FLENS_LAPACK_IMPL_ORGLQ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- orglq ---------------------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orglq_impl(IndexType k, GeMatrix<MA> &A, const DenseVector<VTAU> &tau,
           DenseVector<VWORK> &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    const T Zero(0);

//
//  Perform and apply workspace query
//
    IndexType nb = ilaenv<T>(1, "ORGLQ", "", m, n, k);
    const IndexType lWorkOpt = max(IndexType(1), m) * nb;
    if (work.length()==0) {
        work.resize(max(lWorkOpt,IndexType(1)));
        work(1)=lWorkOpt;
    }
//
//  Quick return if possible
//
    if (m<=0) {
        work(1) = 1;
        return;
    }

    IndexType nbMin = 2;
    IndexType nx = 0;
    IndexType iws = m;

    if ((nb>1) && (nb<k)) {
//
//      Determine when to cross over from blocked to unblocked code.
//
        nx = max(IndexType(0), IndexType(ilaenv<T>(3, "ORGLQ", "", m, n, k)));
        if (nx<k) {
//
//          Determine if workspace is large enough for blocked code.
//
            IndexType ldWork = m;
            iws = ldWork *nb;
            if (work.length()<iws) {
//
//              Not enough workspace to use optimal NB:  reduce NB and
//              determine the minimum value of NB.
//
                nb = work.length() / ldWork;
                nbMin = max(IndexType(2),
                            IndexType(ilaenv<T>(2, "ORGLQ", "", m, n, k)));
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
//      Set A(kk+1:m,1:kk) to zero.
//
        A(_(kk+1,m),_(1,kk)) = Zero;
    } else {
        kk = 0;
    }

//
//  Use unblocked code for the last or only block.
//
    if (kk<n) {
        orgl2(k-kk, A(_(kk+1,m),_(kk+1,n)), tau(_(kk+1, k)), work(_(1,m-kk)));
    }

    if (kk>0) {
        typename GeMatrix<MA>::View Work(m, nb, work);
//
//      Use blocked code
//
        for (IndexType i=ki+1; i>=1; i-=nb) {
            const IndexType ib = min(nb, k-i+1);
            if (i+ib<=m) {
//
//              Form the triangular factor of the block reflector
//              H = H(i) H(i+1) . . . H(i+ib-1)
//
                auto Tr = Work(_(1,ib),_(1,ib)).upper();
                larft(Forward, RowWise,
                      n-i+1,
                      A(_(i,i+ib-1),_(i,n)),
                      tau(_(i,i+ib-1)),
                      Tr);
//
//              Apply H to A(i:m,i+ib:n) from the left
//
                larfb(Right, Trans, Forward, RowWise,
                      A(_(i,i+ib-1),_(i,n)),
                      Tr,
                      A(_(i+ib,m),_(i,n)),
                      Work(_(ib+1,m),_(1,ib)));
            }
//
//          Apply H to rows i:m of current block
//
            orgl2(ib, A(_(i,i+ib-1),_(i,n)), tau(_(i,i+ib-1)), work(_(1,ib)));
//
//          Set columns 1:i-1 of current block to zero
//
            A(_(i,i+ib-1),_(1,i-1)) = Zero;
        }
    }
    work(1) = iws;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- orglq ---------------------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
orglq_impl(IndexType                k,
           GeMatrix<MA>             &A,
           const DenseVector<VTAU>  &tau,
           DenseVector<VWORK>       &work)
{
    cxxlapack::orglq<IndexType>(A.numRows(),
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

//-- orglq ---------------------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VTAU>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
orglq(IndexType k, MA &&A, const VTAU &tau, VWORK &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   ElementType;
    typedef typename RemoveRef<VWORK>::Type VectorWork;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));
    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(tau.length()==k);
    ASSERT((work.length()==0) || (work.length()>=A.numRows()));

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(n>=m);
    ASSERT(m>=k);
    ASSERT(k>=0);
#   endif 

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView    A_org      = A;
    typename VectorWork::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::orglq_impl(k, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename MatrixA::NoView        A_generic      = A;
    typename VectorWork::NoView     work_generic   = work;

    A = A_org;

    if (work_org.length()!=0) {
        work = work_org;
    } else {
        work = ElementType(0);
    }
//
//  Compare results
//
    external::orglq_impl(k, A, tau, work);

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
        std::cerr << "error in: orglq.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: orglq.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ORGLQ_TCC
