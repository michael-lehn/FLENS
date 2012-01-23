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
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LQ_ORML2_TCC
#define FLENS_LAPACK_LQ_ORML2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orml2_generic(Side side, Transpose trans, GeMatrix<MA> &A,
              const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
              DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::IndexType    IndexType;
    typedef typename GeMatrix<MC>::ElementType  T;

    typedef Range<IndexType>    Range;
    const Underscore<IndexType> _;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    const T  One(1);

    const bool noTrans = ((trans==Trans) || (trans==ConjTrans)) ? false
                                                                : true;
//
//  nq is the order of Q
//
    const IndexType nq = (side==Left) ? m : n;

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (k==0)) {
        return;
    }

    IndexType iBeg, iEnd, iInc;
    if ((side==Left && noTrans) || (side==Right && !noTrans)) {
        iBeg = 1;
        iEnd = k;
        iInc = 1;
    } else {
        iBeg = k;
        iEnd = 1;
        iInc = -1;
    }
    iEnd += iInc;

    Range rows = _(1,m);
    Range cols = _(1,n);

    for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
        if (side==Left) {
//
//          H(i) is applied to C(i:m,1:n)
//
            rows = _(i,m);
        } else {
//
//          H(i) is applied to C(1:m,i:n)
//
            cols = _(i,n);
        }
//
//      Apply H(i)
//
        const T Aii = A(i,i);
        A(i,i) = One;
        larf(side, A(i,_(i,nq)), tau(i), C(rows, cols), work);
        A(i,i) = Aii;
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orml2_native(Side side, Transpose trans, GeMatrix<MA> &A,
             const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
             DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::ElementType    T;

    const char      SIDE  = char(side);
    const char      TRANS = getF77LapackChar(trans);
    const INTEGER   M = C.numRows();
    const INTEGER   N = C.numCols();
    const INTEGER   K = A.numRows();
    const INTEGER   LDA = A.leadingDimension();
    const INTEGER   LDC = C.leadingDimension();
    INTEGER         INFO;

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_DECL(dorml2)(&SIDE,
                            &TRANS,
                            &M,
                            &N,
                            &K,
                            A.data(),
                            &LDA,
                            tau.data(),
                            C.data(),
                            &LDC,
                            work.data(),
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO==0);
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orml2(Side side, Transpose trans, GeMatrix<MA> &A,
      const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
      DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::IndexType    IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    const IndexType nq = (side==Left) ? m : n;

    ASSERT(tau.length()==k);

    ASSERT(k<=nq);
    ASSERT(A.numRows()==k);
    ASSERT(A.numCols()==nq);

    if (work.length()>0) {
        if (side==Left) {
            ASSERT(work.length()==n);
        } else {
            ASSERT(work.length()==m);
        }
    }
#   endif
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename GeMatrix<MC>::NoView       C_org      = C;
    typename DenseVector<VWORK>::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    orml2_generic(side, trans, A, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic       = A;
    typename GeMatrix<MC>::NoView       C_generic       = C;
    typename DenseVector<VWORK>::NoView work_generic    = work;

//
//  restore output arguments
//
    A = A_org;
    C = C_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    orml2_native(side, trans, A, tau, C, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: orml2.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: orml2.tcc" << std::endl;
    }
#   endif

}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orml2(Side side, Transpose trans, MA &&A, const VTAU &tau, MC &&C,
      VWORK &&work)
{
    CHECKPOINT_ENTER;
    orml2(side, trans, A, tau, C, work);
    CHECKPOINT_LEAVE;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_LQ_ORML2_TCC
