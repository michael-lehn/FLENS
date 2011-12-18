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
       SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_ORGHR_TCC
#define FLENS_LAPACK_EIG_ORGHR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename IndexType, typename  MA, typename  VTAU>
IndexType
orghr_generic_wsq(IndexType                 iLo,
                  IndexType                 iHi,
                  const GeMatrix<MA>        &A,
                  const DenseVector<VTAU>   &)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const IndexType nh = iHi - iLo;
    const IndexType nb = ilaenv<T>(1, "ORGQR", "", nh, nh, nh);

    return max(IndexType(1), nh)*nb;
}

template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
orghr_generic(IndexType                 iLo,
              IndexType                 iHi,
              GeMatrix<MA>              &A,
              const DenseVector<VTAU>   &tau,
              DenseVector<VW>           &work)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const IndexType nh = iHi - iLo;
    const IndexType nb = ilaenv<T>(1, "ORGQR", "", nh, nh, nh);
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
            A(i,j) = T(0);
        }
        for (IndexType i=j+1; i<=iHi; ++i) {
            A(i,j) = A(i,j-1);
        }
        for (IndexType i=iHi+1; i<=n; ++i) {
            A(i,j) = T(0);
        }
    }
    for (IndexType j=1; j<=iLo; ++j) {
        for (IndexType i=1; i<=n; ++i) {
            A(i,j) = T(0);
        }
        A(j,j) = 1;
    }
    for (IndexType j=iHi+1; j<=n; ++j) {
        for (IndexType i=1; i<=n; ++i) {
            A(i,j) = T(0);
        }
        A(j,j) = T(1);
    }

    if (nh>0) {
//
//      Generate Q(ilo+1:ihi,ilo+1:ihi)
//
        orgqr(nh, A(_(iLo+1,iHi),_(iLo+1,iHi)), tau(_(iLo,iHi-1)), work);
    }
    work(1) = lWorkOpt;
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename IndexType, typename  MA, typename  VTAU>
IndexType
orghr_native_wsq(IndexType                 iLo,
                 IndexType                 iHi,
                 const GeMatrix<MA>        &A,
                 const DenseVector<VTAU>   &tau)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const INTEGER    N      = A.numRows();
    const INTEGER    ILO    = iLo;
    const INTEGER    IHI    = iHi;
    const INTEGER    LDA    = A.leadingDimension();
    T                WORK;
    T                DUMMY;
    const INTEGER    LWORK  = -1;
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dorghr)(&N,
                            &ILO,
                            &IHI,
                            &DUMMY,
                            &LDA,
                            tau.data(),
                            &WORK,
                            &LWORK,
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO==0);
    return WORK;
}


template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
orghr_native(IndexType                 iLo,
             IndexType                 iHi,
             GeMatrix<MA>              &A,
             const DenseVector<VTAU>   &tau,
             DenseVector<VW>           &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const INTEGER    N      = A.numRows();
    const INTEGER    ILO    = iLo;
    const INTEGER    IHI    = iHi;
    const INTEGER    LDA    = A.leadingDimension();
    const INTEGER    LWORK  = work.length();
    INTEGER          INFO;

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dorghr)(&N,
                            &ILO,
                            &IHI,
                            A.data(),
                            &LDA,
                            tau.data(),
                            work.data(),
                            &LWORK,
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO==0);
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename  MA, typename  VTAU>
IndexType
orghr_wsq(IndexType                 iLo,
          IndexType                 iHi,
          const GeMatrix<MA>        &A,
          const DenseVector<VTAU>   &tau)
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
    IndexType ws = orghr_generic_wsq(iLo, iHi, A, tau);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _ws = orghr_native_wsq(iLo, iHi, A, tau);

    if (ws!=_ws) {
        std::cerr << "CXXLAPACK:  ws = " << ws << std::endl;
        std::cerr << "F77LAPACK: _ws = " << _ws << std::endl;
        ASSERT(0);
    }
#   endif

    return ws;
}

template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
orghr(IndexType                 iLo,
      IndexType                 iHi,
      GeMatrix<MA>              &A,
      const DenseVector<VTAU>   &tau,
      DenseVector<VW>           &work)
{
    LAPACK_DEBUG_OUT("orghr");

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
    typename GeMatrix<MA>::NoView       _A    = A;
    typename DenseVector<VW>::NoView    _work = work;
#   endif

//
//  Call implementation
//
    orghr_generic(iLo, iHi, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    if (_work.length()==0) {
        _work.resize(work.length());
    }
    orghr_native(iLo, iHi, _A, tau, _work);

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
        std::cerr << "error in: hrd.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: hrd.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename  MA, typename  VTAU>
IndexType
orghr_wsq(IndexType     iLo,
          IndexType     iHi,
          const MA      &&A,
          const VTAU    &tau)
{
    return orghr_wsq(iLo, iHi, A, tau);
}

template <typename IndexType, typename  MA, typename  VTAU, typename VW>
void
orghr(IndexType     iLo,
      IndexType     iHi,
      MA            &&A,
      const VTAU    &tau,
      VW            &&work)
{
    orghr(iLo, iHi, A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_ORGHR_TCC
