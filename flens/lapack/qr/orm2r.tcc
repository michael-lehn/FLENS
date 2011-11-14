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
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_QR_ORM2R_TCC
#define FLENS_LAPACK_QR_ORM2R_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orm2r(Side side, Transpose trans, MA &&A, const VTAU &tau, MC &&C,
      VWORK &&work)
{
    orm2r(side, trans, A, tau, C, work);
}

//-- ormqr ---------------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
orm2r(Side side, Transpose trans, GeMatrix<MA> &A,
      const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
      DenseVector<VWORK> &work)
{
#   ifndef NDEBUG
    if ((side==Left) && (work.length()<C.numCols())) {
        ASSERT(0);
    }
    if ((side==Right) && (work.length()<C.numRows())) {
        ASSERT(0);
    }
#   endif

    typedef typename GeMatrix<MC>::IndexType    IndexType;
    typedef typename GeMatrix<MC>::ElementType  T;

    typedef Range<IndexType>    Range;
    const Underscore<IndexType> _;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numCols();

    const bool noTrans = ((trans==Trans) || (trans==ConjTrans)) ? false
                                                                : true;
//
//  nq is the order of Q
//
    const IndexType nq = (side==Left) ? m : n;
    
    ASSERT(A.numRows()==nq);
    ASSERT(k<=nq);

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (k==0)) {
        return;
    }

    IndexType iBeg, iEnd, iInc;
    if (((side==Left) && !noTrans) || ((side==Right) && noTrans))
    {
        iBeg = 1;
        iInc = 1;
        iEnd = k+iInc;
    } else {
        iBeg = k;
        iInc = -1;
        iEnd = 1+iInc;
    }

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
        A(i,i) = T(1);
        larf(side, A(_(i,nq), i), tau(i), C(rows, cols), work);
        A(i,i) = Aii;
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_ORM2R_TCC
