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
      SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_EIG_HD2_TCC
#define FLENS_LAPACK_EIG_HD2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename IndexType, typename MA, typename VTAU, typename VW>
void
hd2_generic(IndexType           iLo,
            IndexType           iHi,
            GeMatrix<MA>        &A,
            DenseVector<VTAU>   &tau,
            DenseVector<VW>     &work)
{
    using lapack::larf;
    using lapack::larfg;
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType n = A.numCols();

    for (IndexType i=iLo; i<iHi; ++i) {
//
//      Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
//
        larfg(iHi-i, A(i+1,i), A(_(min(i+2,n),iHi),i), tau(i));

        const T Aii = A(i+1,i);
        A(i+1,i) = 1; 
//
//      Apply H(i) to A(1:ihi,i+1:ihi) from the right
//
        larf(Right, A(_(i+1,iHi),i), tau(i), A(_(1,iHi),_(i+1,iHi)), work);
//
//      Apply H(i) to A(i+1:ihi,i+1:n) from the left
//
        larf(Left, A(_(i+1,iHi),i), tau(i), A(_(i+1,iHi),_(i+1,n)), work);

        A(i+1,i) = Aii;
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
hd2_native(IndexType iLo, IndexType iHi, GeMatrix<MA> &A,
           DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    // TODO: add this assertion to other ..._native implementations
    assert(A.order()==ColMajor);
    typedef typename  GeMatrix<MA>::ElementType     T;

    const INTEGER N     = A.numCols();
    const INTEGER ILO   = iLo;
    const INTEGER IHI   = iHi;
    const INTEGER LDA   = A.leadingDimension();

    INTEGER INFO;

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_IMPL(dgehd2)(&N,
                            &ILO,
                            &IHI,
                            A.data(),
                            &LDA,
                            tau.data(),
                            work.data(),
                            &INFO);
    } else {
        ASSERT(0);
    }

    ASSERT(INFO==0);
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename MA, typename VTAU, typename VW>
void
hd2(IndexType           iLo,
    IndexType           iHi,
    GeMatrix<MA>        &A,
    DenseVector<VTAU>   &tau,
    DenseVector<VW>     &work)
{
    LAPACK_DEBUG_OUT("hd2");

    using std::max;
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
    ASSERT(tau.firstIndex()<=iLo);
    ASSERT(tau.lastIndex()>=iHi-1);
    ASSERT(tau.inc()>0);
    ASSERT(work.length()>=A.numRows());

    ASSERT(1<=iLo);
    ASSERT(iLo<=iHi);
    ASSERT(iHi<=max(IndexType(1), A.numCols()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       _A = A;
    typename DenseVector<VTAU>::NoView  _tau = tau;
    typename DenseVector<VW>::NoView    _work = work;
#   endif

//
//  Call implementation
//
    hd2_generic(iLo, iHi, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    hd2_native(iLo, iHi, _A, _tau, _work);

    if (! isIdentical(A, _A, "A", "A_")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        ASSERT(0);
    }

    if (! isIdentical(tau, _tau, "tau", "tau_")) {
        std::cerr << "CXXLAPACK:  tau = " << tau << std::endl;
        std::cerr << "F77LAPACK: _tau = " << _tau << std::endl;
        ASSERT(0);
    }

    if (! isIdentical(work, _work, "work", "work_")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MA, typename VTAU, typename VW>
void
hd2(IndexType   iLo,
    IndexType   iHi,
    MA          &&A,
    VTAU        &&tau,
    VW          &&work)
{
    hd2(iLo, iHi, A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_HD2_TCC
