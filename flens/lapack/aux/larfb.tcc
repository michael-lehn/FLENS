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
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_AUX_LARFB_TCC
#define FLENS_LAPACK_AUX_LARFB_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb(Side side, Transpose transH, Direction direction, StoreVectors storeV,
      const MV &V, const MT &T, MC &&C, MWORK &&Work)
{
    larfb(side, transH, direction, storeV, V, T, C, Work);
}

//-- larf ----------------------------------------------------------------------
template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb(Side side, Transpose transH, Direction direction, StoreVectors storeV,
      const GeMatrix<MV> &V, const TrMatrix<MT> &T, GeMatrix<MC> &C,
      GeMatrix<MWORK> &Work)
{
    ASSERT(transH!=Conj);
    ASSERT(Work.numRows()>= (side==Left) ? C.numCols() : C.numRows());
    ASSERT(Work.numCols()>=T.dim());

    typedef typename GeMatrix<MC>::ElementType  ElementType;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const ElementType   One(1);

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = T.dim();

//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return;
    }

    const Transpose transT = (transH==NoTrans) ? ConjTrans : NoTrans;

    if (storeV==ColumnWise) {
        if (direction==Forward) {
//
//          Let  V =  ( V1 )    (first K rows)
//                    ( V2 )
//          where  V1  is unit lower triangular.
//
            if (side==Left) {
//
//              Form  H * C  or  H**T * C  where  C = ( C1 )
//                                                    ( C2 )
//
                IndexType lastV = max(k, ilalr(V));
                const auto V1 = V(_(1,k),_);
                const auto V2 = V(_(k+1,lastV),_);

                IndexType lastC = ilalc(C(_(1,lastV),_));
                auto C1 = C(_(1,k),_(1,lastC));
                auto C2 = C(_(k+1,lastV),_(1,lastC));

//
//              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
//
//              W := C1**T
//
                //std::cerr << "-> Work = " << Work << std::endl;
                auto W = Work(_(1,lastC),_(1,k));

                blas::copy(Trans, C1, W);
                //std::cerr << "-> Work = " << Work << std::endl;
                //std::cerr << "-> V = " << V << std::endl;
                //std::cerr << "-> C = " << C << std::endl;
                //std::cerr << "-> C1 = " << C1 << std::endl;
                //std::cerr << "-> W    = " << W << std::endl;
//
//              W := W * V1
//
                blas::mm(Right, NoTrans, One, V1.lowerUnit(), W);

                if (lastV>k) {
//
//                  W := W + C2**T *V2
//
                    blas::mm(Trans, NoTrans, One, C2, V2, One, W);
                }
//
//              W := W * T**T  or  W * T
//
                blas::mm(Right, transT, One, T, W);
//
//              C := C - V * W**T
//
                if (lastV>k) {
//
//                  C2 := C2 - V2 * W**T
//
                    blas::mm(NoTrans, Trans, -One, V2, W, One, C2);
                }
//
//              W := W * V1**T
//
                blas::mm(Right, Trans, One, V1.lowerUnit(), W);
//
//              C1 := C1 - W**T
//
                blas::axpy(Trans, -One, W, C1);
            } else if (side==Right) {
//
//              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
//
                IndexType lastV = max(k, ilalr(V));
                const auto V1 = V(_(1,k),_);
                const auto V2 = V(_(k+1,lastV),_);

                IndexType lastC = ilalr(C(_,_(1,lastV)));
                auto C1 = C(_(1,lastC),_(1,k));
                auto C2 = C(_(1,lastC),_(k+1,lastV));
//
//              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
//
//              W := C1
//
                auto W = Work(_(1,lastC),_(1,k));
                blas::copy(NoTrans, C1, W);
//
//              W := W * V1
//
                blas::mm(Right, NoTrans, One, V1.lowerUnit(), W);

                if (lastV>k) {
//
//                  W := W + C2 * V2
//
                    blas::mm(NoTrans, NoTrans, One, C2, V2, One, W);
                }
//
//              W := W * T  or  W * T**T
//
                blas::mm(Right, transH, One, T, W);
//
//              C := C - W * V**T
//
                if (lastV>k) {
//
//                  C2 := C2 - W * V2**T
//
                    blas::mm(NoTrans, Trans, -One, W, V2, One, C2);
                }
//
//              W := W * V1**T
//
                blas::mm(Right, Trans, One, V1.lowerUnit(), W);
//
//              C1 := C1 - W
//
                blas::axpy(NoTrans, -One, W, C1);

                //std::cerr << "8) C1 = " << C1 << std::endl;
            }
        } else if (direction==Backward) {
            // Lehn: I will implement it as soon as someone needs it
            ASSERT(0);
        }
    } else if (storeV==RowWise) {
        // Lehn: I will implement it as soon as someone needs it
        ASSERT(0);
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFB_TCC
