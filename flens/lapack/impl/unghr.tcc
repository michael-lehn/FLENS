/*
 *   Copyright (c) 2013, Michael Lehn
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
       SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_UNGHR_TCC
#define FLENS_LAPACK_IMPL_UNGHR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename  MA, typename  VTAU>
IndexType
unghr_wsq_impl(IndexType                 iLo,
               IndexType                 iHi,
               const GeMatrix<MA>        &A,
               const DenseVector<VTAU>   &)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const IndexType nh = iHi - iLo;
    const IndexType nb = ilaenv<T>(1, "UNGQR", "", nh, nh, nh);

    return max(IndexType(1), nh)*nb;
}

template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
unghr_impl(IndexType                 iLo,
           IndexType                 iHi,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VW>           &work)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;

    const T Zero(0), One(1);

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const IndexType nh = iHi - iLo;
    const IndexType nb = ilaenv<T>(1, "UNGQR", "", nh, nh, nh);
    const IndexType lWorkOpt = max(IndexType(1), nh)*nb;

//
//  Apply workspace query if needed
//
    if (work.length()==0) {
        work.resize(lWorkOpt);
    }
//
//  Quick return if possible
//
    if (n==0) {
        work(1) = 1;
    }
//
//  Shift the vectors which define the elementary reflectors one
//  column to the right, and set the first ilo and the last n-ihi
//  rows and columns to those of the unit matrix
//
    for (IndexType j=iHi; j>iLo; --j) {
        for (IndexType i=1; i<=j-1; ++i) {
            A(i,j) = Zero;
        }
        for (IndexType i=j+1; i<=iHi; ++i) {
            A(i,j) = A(i,j-1);
        }
        for (IndexType i=iHi+1; i<=n; ++i) {
            A(i,j) = Zero;
        }
    }
    for (IndexType j=1; j<=iLo; ++j) {
        for (IndexType i=1; i<=n; ++i) {
            A(i,j) = Zero;
        }
        A(j,j) = One;
    }
    for (IndexType j=iHi+1; j<=n; ++j) {
        for (IndexType i=1; i<=n; ++i) {
            A(i,j) = Zero;
        }
        A(j,j) = One;
    }

    if (nh>0) {
//
//      Generate Q(ilo+1:ihi,ilo+1:ihi)
//
        ungqr(A(_(iLo+1,iHi),_(iLo+1,iHi)), tau(_(iLo,iHi-1)), work);
    }
    work(1) = lWorkOpt;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename  MA, typename  VTAU>
IndexType
unghr_wsq_impl(IndexType                 iLo,
               IndexType                 iHi,
               const GeMatrix<MA>        &A,
               const DenseVector<VTAU>   &tau)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;

    ElementType      WORK;
    ElementType      DUMMY;
    const IndexType  LWORK  = -1;

    cxxlapack::unghr<IndexType>(A.numRows(),
                                iLo,
                                iHi,
                                &DUMMY,
                                A.leadingDimension(),
                                tau.data(),
                                &WORK,
                                LWORK);
    return WORK;
}


template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
unghr_impl(IndexType                 iLo,
           IndexType                 iHi,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VW>           &work)
{
    cxxlapack::unghr<IndexType>(A.numRows(),
                                iLo,
                                iHi,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename  MA, typename  VTAU>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value,
         IndexType>::Type
unghr_wsq(IndexType                 iLo,
          IndexType                 iHi,
          const MA                  &A,
          const VTAU                &tau)
{
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numCols();
    if (n==0) {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    } else {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    }
    ASSERT(tau.length()==(n-1));
#   endif

//
//  Call implementation
//
    IndexType ws = LAPACK_SELECT::unghr_wsq_impl(iLo, iHi, A, tau);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _ws = external::unghr_wsq_impl(iLo, iHi, A, tau);

    if (ws!=_ws) {
        std::cerr << "CXXLAPACK:  ws = " << ws << std::endl;
        std::cerr << "F77LAPACK: _ws = " << _ws << std::endl;
        ASSERT(0);
    }
#   endif

    return ws;
}

template <typename IndexType, typename  MA, typename  VTAU, typename VW>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VW>::value,
         void>::Type
unghr(IndexType                     iLo,
      IndexType                     iHi,
      MA                            &&A,
      const VTAU                    &tau,
      VW                            &&work)
{
    LAPACK_DEBUG_OUT("unghr");
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VW>::Type    VectorW;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numCols();
    if (n==0) {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    } else {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    }
    ASSERT(tau.length()==(n-1));
    ASSERT(work.length()>=iHi-iLo);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView    _A    = A;
    typename VectorW::NoView    _work = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::unghr_impl(iLo, iHi, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }
    external::unghr_impl(iLo, iHi, _A, tau, _work);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "A_")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: unghr.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: unghr.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNGHR_TCC
