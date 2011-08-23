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
       SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
      $                   LDC, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_ORMHR_TCC
#define FLENS_LAPACK_EIG_ORMHR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

using std::max;

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
void
ormhr(Side side, Transpose trans, IndexType iLo, IndexType iHi,
      MA &&A, const VTAU &tau, MC &&C, VWORK &&work)
{
    ormhr(side, trans, iLo, iHi, A, tau, C, work);
}

//-- ormhr ---------------------------------------------------------------------
template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
void
ormhr_impl(Side side, Transpose trans, IndexType iLo, IndexType iHi,
           GeMatrix<MA> &A, const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
           DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::ElementType  T;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

//
//  nq is the order of Q and nw is the minimum dimension of WORK
//
    const IndexType nq = (side==Left) ? m : n;
    const IndexType nw = (side==Left) ? n : m;

//  TODO: implement a better way for setting up the opt string
    char opt[3];
    opt[0] = (Side==Left) ? 'L' : 'R';
    if (trans==NoTrans) {
        opt[1] = 'N';
    } else if (trans==Trans) {
        opt[1] = 'T';
    } else if (trans==ConjTrans) {
        opt[1] = 'C';
    } else if (trans==Conj) {
        opt[1] = 'R';
    }

    IndexType nb;
    if (side==Left) {
        nb = ilaenv<T>(1, "ORMQR", opt, nh, n, nh);
    } else {
        nb = ilaenv<T>(1, "ORMQR", opt, m, nh, nh);
    }
    IndexType lWorkOpt = max(IndexType(1), nw)*nb;
    work(1) = lWorkOpt;

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (nh==0)) {
        work(1) = 1;
        return;
    }

    if (side==Left) {
        mi = nh;
        ni = n;
        i1 = iLo + 1;
        i2 = 1;
    } else {
        mi = m
        ni = nh
        i1 = 1
        i2 = ilo + 1
    }

    auto _A     = A(_(iLo+1,iLo+mi),_(iLo,iLo+nh-1));
    auto _tau   = tau(_(iLo,iLo+nh-1));
    auto _C     = C(_(i1,i1+mi-1),_(i2,i2+ni-1));

    ormqr(side, trans, A(_(iLo+1,iLo+mi),_(iLo,iLo+ni-1)), tau(), C(), work)

    work(1) = lWorkOpt;
}

template <typename IndexType, typename  MA, typename  VTAU,
          typename  MC, typename  VWORK>
void
ormhr(Side side, Transpose trans, IndexType iLo, IndexType iHi,
      GeMatrix<MA> &A, const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
      DenseVector<VWORK> &work)
{
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(m==1);
    ASSERT(n==1);
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

    ormhrImpl(side, trans, iLo, iHi, A, tau, C, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_ORMHR_TCC
