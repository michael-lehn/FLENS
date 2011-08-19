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
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
      IMPLICIT NONE
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_AUX_LARFT_TCC
#define FLENS_LAPACK_AUX_LARFT_TCC 1

#include <flens/blas/blas.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MV, typename VTAU, typename MT>
void
larft(Direction direction, StoreVectors storeVectors, IndexType n,
      MV &&V, const VTAU &tau, MT &&T)
{
    larft(direction, storeVectors, n, V, tau, T);
}

//-- larf ----------------------------------------------------------------------
template <typename IndexType, typename MV, typename VTAU, typename MT>
void
larft(Direction direction, StoreVectors storeVectors, IndexType n,
      GeMatrix<MV> &V, const DenseVector<VTAU> &tau, TrMatrix<MT> &T)
{
    typedef typename GeMatrix<MV>::ElementType  ElementType;
    const Underscore<IndexType> _;

//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    const IndexType k = T.dim();

//  Lehn: as long as we do not have col-views for TrMatrix we get them
//        via a GeMatrixView
    auto _T = T.general();

    if (direction==Forward) {
        IndexType lastV;
        IndexType prevLastV = n;
        for (IndexType i=1; i<=k; ++i) {
            prevLastV = max(i, prevLastV);
            if (tau(i)==ElementType(0)) {
//
//              H(i)  =  I
//
                for (IndexType j=1; j<=i; ++j) {
                    T(j,i) = 0;
                }
            } else {
//
//              general case
//
                ElementType Vii = V(i,i);
                V(i,i) = 1;
                if (storeVectors==ColumnWise) {
//                  Skip any trailing zeros.
                    for (lastV=n; lastV>=i+1; --lastV) {
                        if (V(lastV,i)!=ElementType(0)) {
                            break;
                        }
                    }
                    IndexType j = min(lastV, prevLastV);
//
//                  T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
//
                    blas::mv(Trans, -tau(i),
                             V(_(i,j),_(1,i-1)), V(_(i,j),i),
                             ElementType(0),
                             _T(_(1,i-1),i));
                } else { /* storeVectors==RowWise */
                    // Lehn: I will implement it as soon as someone needs it
                    ASSERT(0);
                }
                V(i,i) = Vii;
//
//              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
//
                blas::mv(NoTrans,
                         _T(_(1,i-1),_(1,i-1)).upper(),
                         _T(_(1,i-1),i));
                T(i,i) = tau(i);
                if (i>1) {
                    prevLastV = max(prevLastV, lastV);
                } else {
                    prevLastV = lastV;
                }
            }
        }
    } else {
        // Lehn: I will implement it as soon as someone needs this case
        ASSERT(0);
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFT_TCC
