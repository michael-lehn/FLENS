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
      SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_UNMHR_TCC
#define FLENS_LAPACK_IMPL_UNMHR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename  MC>
IndexType
unmhr_wsq_impl(Side                      side,
               Transpose                 trans,
               IndexType                 iLo,
               IndexType                 iHi,
               const GeMatrix<MC>        &C)
{
    using std::max;

    typedef typename GeMatrix<MC>::ElementType  T;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType nh = iHi - iLo;

//
//  nw is the minimum dimension of WORK
//
    const IndexType nw = (side==Left) ? n : m;

//  TODO: implement a better way for setting up the opt string
    char opt[3];
    opt[0] = getF77Char(side);
    opt[1] = getF77Char(trans);
    opt[2] = 0;

    IndexType nb;
    if (side==Left) {
        nb = ilaenv<T>(1, "UNMQR", opt, nh, n, nh);
    } else {
        nb = ilaenv<T>(1, "UNMQR", opt, m, nh, nh);
    }
    return max(IndexType(1), nw)*nb;
}

template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
void
unmhr_impl(Side                      side,
           Transpose                 trans,
           IndexType                 iLo,
           IndexType                 iHi,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           GeMatrix<MC>              &C,
           DenseVector<VWORK>        &work)
{
    using std::max;

    typedef typename GeMatrix<MC>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType nh = iHi - iLo;

//
//  nw is the minimum dimension of WORK
//
    const IndexType nw = (side==Left) ? n : m;

//  TODO: implement a better way for setting up the opt string
    char opt[3];
    opt[0] = getF77Char(side);
    opt[1] = getF77Char(trans);
    opt[2] = 0;

    IndexType nb;
    if (side==Left) {
        nb = ilaenv<T>(1, "UNMQR", opt, nh, n, nh);
    } else {
        nb = ilaenv<T>(1, "UNMQR", opt, m, nh, nh);
    }
    IndexType lWorkOpt = max(IndexType(1), nw)*nb;

//
//  Apply worksize query
//
    if (work.length()==0) {
        work.resize(lWorkOpt);
    }

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (nh==0)) {
        work(1) = 1;
        return;
    }

    const auto tau_ = tau(_(iLo,iHi-1));
    auto       A_  = A(_(iLo+1,iHi),_(iLo,iHi-1));

    if (side==Left) {

        auto C_  = C(_(iLo+1,iHi),_);

        unmqr(Left, trans, A_, tau_, C_, work);

    } else {

        auto C_  = C(_,_(iLo+1,iHi));

        unmqr(Right, trans, A_, tau_, C_, work);
    }

    work(1) = lWorkOpt;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename  MC>
IndexType
unmhr_wsq_impl(Side                      side,
               Transpose                 trans,
               IndexType                 iLo,
               IndexType                 iHi,
               const GeMatrix<MC>        &C)
{
    using std::max;

    typedef typename GeMatrix<MC>::ElementType  ElementType;

    IndexType m = C.numRows();
    IndexType n = C.numCols();

    const IndexType     LDA    = (side==Left) ? max(IndexType(1), m)
                                              : max(IndexType(1), n);
    ElementType         WORK, DUMMY;
    const IndexType     LWORK  = -1;

    cxxlapack::unmhr<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                iLo,
                                iHi,
                                &DUMMY,
                                LDA,
                                &DUMMY,
                                &DUMMY,
                                C.leadingDimension(),
                                &WORK,
                                LWORK);
    return real(WORK);
}

template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
void
unmhr_impl(Side                      side,
           Transpose                 trans,
           IndexType                 iLo,
           IndexType                 iHi,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           GeMatrix<MC>              &C,
           DenseVector<VWORK>        &work)
{
    cxxlapack::unmhr<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                iLo,
                                iHi,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                C.data(),
                                C.leadingDimension(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename  MC>
IndexType
unmhr_wsq(Side                      side,
          Transpose                 trans,
          IndexType                 iLo,
          IndexType                 iHi,
          const GeMatrix<MC>        &C)
{
    LAPACK_DEBUG_OUT("unmhr_wsq");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    if (side==Left) {
        if (m==0) {
            ASSERT(iLo==1);
            ASSERT(iHi==0);
        } else {
            ASSERT(1<=iLo);
            ASSERT(iLo<=iHi);
            ASSERT(iHi<=m);
        }
    }
    if (side==Right) {
        if (n==0) {
            ASSERT(iLo==1);
            ASSERT(iHi==0);
        } else {
            ASSERT(1<=iLo);
            ASSERT(iLo<=iHi);
            ASSERT(iHi<=n);
        }
    }
#   endif

//
//  Call implementation
//
    IndexType ws = LAPACK_SELECT::unmhr_wsq_impl(side, trans, iLo, iHi, C);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    IndexType ws_ = external::unmhr_wsq_impl(side, trans, iLo, iHi, C);

    if (ws!=ws_) {
        std::cerr << "CXXLAPACK:  ws = " << ws << std::endl;
        std::cerr << "F77LAPACK: ws_ = " << ws_ << std::endl;
        ASSERT(0);
    }
#   endif
    return ws;
}

template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexGeMatrix<MC>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
unmhr(Side                          side,
      Transpose                     trans,
      IndexType                     iLo,
      IndexType                     iHi,
      MA                            &&A,
      const VTAU                    &&tau,
      MC                            &&C,
      VWORK                         &&work)
{
    LAPACK_DEBUG_OUT("unmhr");

    using std::max;

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<MC>::Type        MatrixC;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    if (side==Left) {
        if (m==0) {
            ASSERT(iLo==1);
            ASSERT(iHi==0);
        } else {
            ASSERT(1<=iLo);
            ASSERT(iLo<=iHi);
            ASSERT(iHi<=m);
        }
    }
    if (side==Right) {
        if (n==0) {
            ASSERT(iLo==1);
            ASSERT(iHi==0);
        } else {
            ASSERT(1<=iLo);
            ASSERT(iLo<=iHi);
            ASSERT(iHi<=n);
        }
    }
    ASSERT((side==Left  && A.numCols()==m)
        || (side==Right && A.numCols()==n));
    ASSERT((side==Left  && tau.length()==(m-1))
        || (side==Right && tau.length()==(n-1)));
    ASSERT((side==Left  && work.length()>=max(IndexType(1),n))
        || (side==Right && work.length()>=max(IndexType(1),m)));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView       A_org      = A;
    typename MatrixC::NoView       C_org      = C;
    typename VectorWork::NoView    work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::unmhr_impl(side, trans, iLo, iHi, A, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView       A_generic       = A;
    typename MatrixC::NoView       C_generic       = C;
    typename VectorWork::NoView    work_generic    = work;

//
//  restore output arguments
//
    A = A_org;
    C = C_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::unmhr_impl(side, trans, iLo, iHi, A, tau, C, work);

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
        std::cerr << "error in: unmhr.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: unmhr.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNMHR_TCC
