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
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 *
 */

#ifndef FLENS_LAPACK_QR_ORMQR_TCC
#define FLENS_LAPACK_QR_ORMQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormqr(Side side, Transpose &trans, MA &&A, const VTAU &tau, MC &&C,
      VWORK &&work)
{
    ormqr(side, trans, A, tau, C, work);
}

//-- ormqr ---------------------------------------------------------------------
template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormqr(Side side, Transpose &trans, GeMatrix<MA> &A,
      const DenseVector<VTAU> &tau, GeMatrix<MC> &C,
      DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    typedef typename GeMatrix<MC>::View         GeView;
    typedef typename GeView::Engine             GeViewEngine;

    const Underscore<IndexType> _;


//
//  Paramter for maximum block size and buffer for TrMatrix Tr.
//
    const IndexType     nbMax = 64;
    const IndexType     ldt = nbMax + 1;
    T                   trBuffer[nbMax*ldt];

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numCols();

    ASSERT(tau.length()==k);

    const bool noTrans = ((trans==Trans) || (trans==ConjTrans)) ? false
                                                                : true;
//
//  nq is the order of Q and nw is the minimum dimension of work
//
    IndexType nq, nw;
    if (side==Left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
//
//  Determine the block size.  nb may be at most nbMax, where nbMax
//  is used to define the local array tr.
//
    char opt[3];
    opt[0] = (side==Left) ? 'L' : 'R';
    if (trans==NoTrans) {
        opt[1] = 'N';
    } else if (trans==Conj) {
        opt[1] = 'R';
    } else if (trans==Trans) {
        opt[1] = 'T';
    } else if (trans==ConjTrans) {
        opt[1] = 'C';
    }
    opt[2] = 0;

    IndexType nb = min(nbMax, IndexType(ilaenv<T>(1, "ORMQR", opt, m, n, k)));
    IndexType lWorkOpt = max(IndexType(1), nw)*nb;

    if (work.length()==0) {
        work.resize(lWorkOpt);
    }

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (k==0)) {
        work(1) = 1;
        return;
    }

    IndexType nbMin = 2;
    IndexType iws;
    if ((nb>1) && (nb<k)) {
        iws = lWorkOpt;
        if (work.length()<iws) {
            nb = work.length()/nw;
            nbMin = max(nbMin, IndexType(ilaenv<T>(2, "ORMQR", opt, m, n, k)));
        }
    } else {
        iws = nw;
    }

    if ((nb<nbMin) || (nb>=k)) {
//
//      Use unblocked code
//
        orm2r(side, trans, A, tau, C, work);
    } else {
//
//      Use blocked code
//
        IndexType iBeg, iInc, iEnd;
        if (((side==Left) && (!noTrans)) || ((side==Right) && (noTrans))) {
            iBeg = 1;
            iInc = nb;
            iEnd = ((k-1)/nb)*nb + 1 + iInc;
        } else {
            iBeg = ((k-1)/nb)*nb + 1;
            iInc = -nb;
            iEnd = 1 + iInc;
        }

        IndexType ic, jc;
        if (side==Left) {
            jc = 1;
        } else {
            ic = 1;
        }

        typename GeMatrix<MA>::View Work(nw, nb, work);

        for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
            const IndexType ib = min(nb, k-i+1);
            GeView          Tr = GeViewEngine(ib, ib, trBuffer, ib);
            
//
//          Form the triangular factor of the block reflector
//          H = H(i) H(i+1) . . . H(i+ib-1)
//
            larft(Forward, ColumnWise, nq-i+1,
                  A(_(i,nq),_(i,i+ib-1)), tau(_(i,i+ib-1)), Tr.upper());

            if (side==Left) {
//
//              H or H**T is applied to C(i:m,1:n)
//
                ic = i;
            } else {
//
//              H or H**T is applied to C(1:m,i:n)
//
                jc = i;
            }

            larfb(side, trans,
                  Forward, ColumnWise,
                  A(_(i,nq),_(i,i+ib-1)), Tr.upper(),
                  C(_(ic,m),_(jc,n)),
                  Work);
        }
    }
    work(1) = lWorkOpt;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_QR_ORMQR_TCC
