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

//== generic lapack implementation =============================================

template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb_generic(Side                  side,
              Transpose             transH,
              Direction             direction,
              StoreVectors          storeVectors,
              const GeMatrix<MV>    &V,
              const TrMatrix<MT>    &Tr,
              GeMatrix<MC>          &C,
              GeMatrix<MWORK>       &Work)
{
    using lapack::ilalc;
    using lapack::ilalr;
    using std::max;

    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const T                     One(1);

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = Tr.dim();

//
//  Quick return if possible
//
    if ((m==0) || (n==0)) {
        return;
    }

    const Transpose transT = (transH==NoTrans) ? ConjTrans : NoTrans;

    if (storeVectors==ColumnWise) {
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
                const IndexType lastV = max(k, ilalr(V));
                const auto V1 = V(_(1,k),_);
                const auto V2 = V(_(k+1,lastV),_);

                const IndexType lastC = ilalc(C(_(1,lastV),_));
                auto C1 = C(_(1,k),_(1,lastC));
                auto C2 = C(_(k+1,lastV),_(1,lastC));

//
//              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
//
//              W := C1**T
//
                auto W = Work(_(1,lastC),_(1,k));
                blas::copy(Trans, C1, W);
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
                blas::mm(Right, transT, One, Tr, W);
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
                const IndexType lastV = max(k, ilalr(V));
                const auto V1 = V(_(1,k),_);
                const auto V2 = V(_(k+1,lastV),_);

                const IndexType lastC = ilalr(C(_,_(1,lastV)));
                auto C1 = C(_(1,lastC),_(1,k));
                auto C2 = C(_(1,lastC),_(k+1,lastV));
//
//              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
//
//              W := C1
//
                auto W = Work(_(1,lastC),_(1,k));
                W = C1;
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
                blas::mm(Right, transH, One, Tr, W);
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
            }
        } else if (direction==Backward) {
            // Lehn: I will implement it as soon as someone needs it
            ASSERT(0);
        }
    } else if (storeVectors==RowWise) {

        if (direction==Forward) {
//
//          Let  V =  ( V1  V2 )    (V1: first K columns)
//          where  V1  is unit upper triangular.
//
            if (side==Left) {
//
//              Form  H * C  or  H**T * C  where  C = ( C1 )
//                                                    ( C2 )
//
                const IndexType lastV = max(k, ilalc(V));
                const auto V1 = V(_,_(1,k));
                const auto V2 = V(_,_(k+1,lastV));

                const IndexType lastC = ilalc(C(_(1,lastV),_));
                auto C1 = C(_(  1,    k),_(1,lastC));
                auto C2 = C(_(k+1,lastV),_(1,lastC));
//
//              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T)
//                                                              (stored in WORK)
//              W := C1**T
//
                auto W = Work(_(1,lastC),_(1,k));
                blas::copy(Trans, C1, W);
//
//              W := W * V1**T
//
                blas::mm(Right, Trans, One, V1.upperUnit(), W);

                if (lastV>k) {
//
//                  W := W + C2**T*V2**T
//
                    blas::mm(Trans, Trans, One, C2, V2, One, W);
                }
//
//              W := W * T**T  or  W * T
//
                blas::mm(Right, transT, One, Tr, W);
//
//              C := C - V**T * W**T
//
                if (lastV>k) {
//
//                  C2 := C2 - V2**T * W**T
//
                    blas::mm(Trans, Trans, -One, V2, W, One, C2);
                }
//
//              W := W * V1
//
                blas::mm(Right, NoTrans, One, V1.upperUnit(), W);
//
//              C1 := C1 - W**T
//
                blas::axpy(Trans, -One, W, C1);
            } else if (side==Right) {
//
//              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
//
                const IndexType lastV = max(k, ilalc(V));
                const auto V1 = V(_,_(1,k));
                const auto V2 = V(_,_(k+1,lastV));

                const IndexType lastC = ilalr(C(_,_(1,lastV)));
                auto C1 = C(_(1,lastC),_(1,k));
                auto C2 = C(_(1,lastC),_(k+1,lastV));
//
//              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
//
//              W := C1
//
                auto W = Work(_(1,lastC),_(1,k));
                W = C1;
//
//              W := W * V1**T
//
                blas::mm(Right, Trans, One, V1.upperUnit(), W);

                if (lastV>k) {
//
//                  W := W + C2 * V2**T
//
                    blas::mm(NoTrans, Trans, One, C2, V2, One, W);
                }
//
//              W := W * T  or  W * T**T
//
                blas::mm(Right, transH, One, Tr, W);
//
//              C := C - W * V
//
                if (lastV>k) {
//
//                  C2 := C2 - W * V2
//
                    blas::mm(NoTrans, NoTrans, -One, W, V2, One, C2);
                }
//
//              W := W * V1
//
                blas::mm(Right, NoTrans, One, V1.upperUnit(), W);
//
//              C1 := C1 - W
//
                blas::axpy(NoTrans, -One, W, C1);
            }
        } else if (direction==Backward) {
            if (side==Left) {
                // Lehn: I will implement it as soon as someone needs it
                ASSERT(0);
            } else if (side==Right) {
                // Lehn: I will implement it as soon as someone needs it
                ASSERT(0);
            }
        }
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb_native(Side                   side,
             Transpose              transH,
             Direction              direction,
             StoreVectors           storeVectors,
             const GeMatrix<MV>     &V,
             const TrMatrix<MT>     &Tr,
             GeMatrix<MC>           &C,
             GeMatrix<MWORK>        &Work)
{
    typedef typename TrMatrix<MT>::ElementType  T;

    const char      SIDE    = char(side);
    const char      TRANS   = getF77LapackChar(transH);
    const char      DIRECT  = char(direction);
    const char      STOREV  = char(storeVectors);
    const INTEGER   M       = C.numRows();
    const INTEGER   N       = C.numCols();
    const INTEGER   K       = Tr.dim();
    const INTEGER   LDV     = V.leadingDimension();
    const INTEGER   LDT     = Tr.leadingDimension();
    const INTEGER   LDC     = C.leadingDimension();
    const INTEGER   LDWORK  = Work.leadingDimension();

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dlarfb)(&SIDE,
                            &TRANS,
                            &DIRECT,
                            &STOREV,
                            &M,
                            &N,
                            &K,
                            V.data(),
                            &LDV,
                            Tr.data(),
                            &LDT,
                            C.data(),
                            &LDC,
                            Work.data(),
                            &LDWORK);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb(Side                  side,
      Transpose             transH,
      Direction             direction,
      StoreVectors          storeV,
      const GeMatrix<MV>    &V,
      const TrMatrix<MT>    &Tr,
      GeMatrix<MC>          &C,
      GeMatrix<MWORK>       &Work)
{
    LAPACK_DEBUG_OUT("larfb");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(transH!=Conj);

    if (side==Left) {
        ASSERT(Work.numRows()>=C.numCols());
    } else {
        ASSERT(Work.numRows()>=C.numRows());
    }
    ASSERT(Work.numCols()==Tr.dim());
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MC>::NoView       C_org = C;
    typename GeMatrix<MWORK>::NoView    Work_org = Work;
#   endif

//
//  Call implementation
//
    larfb_generic(side, transH, direction, storeV, V, Tr, C, Work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename GeMatrix<MC>::NoView       C_generic = C;
    typename GeMatrix<MWORK>::NoView    Work_generic = Work;

    C    = C_org;
    Work = Work_org;
//
//  Compare results
//
    larfb_native(side, transH, direction, storeV, V, Tr, C, Work);

    bool failed = false;
    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }
    if (! isIdentical(Work_generic, Work, " Work_generic", "Work")) {
        std::cerr << "CXXLAPACK: Work_generic = " << Work_generic << std::endl;
        std::cerr << "F77LAPACK: Work = " << Work << std::endl;
        failed = true;
    }
    if (failed) {
        std::cerr << "side =      " << char(side) << std::endl;
        std::cerr << "transH =    " << transH << std::endl;
        std::cerr << "direction = " << char(direction) << std::endl;
        std::cerr << "storeV =    " << char(storeV) << std::endl;
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MV, typename MT, typename MC, typename MWORK>
void
larfb(Side              side,
      Transpose         transH,
      Direction         direction,
      StoreVectors      storeV,
      const MV          &V,
      const MT          &Tr,
      MC                &&C,
      MWORK             &&Work)
{
    larfb(side, transH, direction, storeV, V, Tr, C, Work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFB_TCC
