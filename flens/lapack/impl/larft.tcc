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

#ifndef FLENS_LAPACK_IMPL_LARFT_TCC
#define FLENS_LAPACK_IMPL_LARFT_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename N, typename MV, typename VTAU, typename MT>
void
larft_impl(Direction direction, StoreVectors storeVectors, N n,
           GeMatrix<MV> &V, const DenseVector<VTAU> &tau, TrMatrix<MT> &Tr)
{
    using std::max;
    using std::min;

    typedef typename TrMatrix<MT>::IndexType    IndexType;
    typedef typename TrMatrix<MT>::ElementType  T;

    const Underscore<IndexType> _;

    const T Zero(0);

//
//  Quick return if possible
//
    if (n==0) {
        return;
    }

    const IndexType k = Tr.dim();

//  Lehn: as long as we do not have col-views for TrMatrix we get them
//        via a GeMatrixView
    auto _Tr = Tr.general();

    if (direction==Forward) {
        IndexType lastV = -1;
        IndexType prevLastV = n;
        for (IndexType i=1; i<=k; ++i) {
            prevLastV = max(i, prevLastV);
            if (tau(i)==Zero) {
//
//              H(i)  =  I
//
                for (IndexType j=1; j<=i; ++j) {
                    Tr(j,i) = Zero;
                }
            } else {
//
//              general case
//
                T Vii = V(i,i);
                V(i,i) = 1;
                if (storeVectors==ColumnWise) {
//                  Skip any trailing zeros.
                    for (lastV=n; lastV>=i+1; --lastV) {
                        if (V(lastV,i)!=Zero) {
                            break;
                        }
                    }
                    IndexType j = min(lastV, prevLastV);
//
//                  T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
//
                    blas::mv(Trans, -tau(i),
                             V(_(i,j),_(1,i-1)), V(_(i,j),i),
                             Zero,
                             _Tr(_(1,i-1),i));
                } else { /* storeVectors==RowWise */
//                  Skip any trailing zeros.
                    for (lastV=n; lastV>=i+1; --lastV) {
                        if (V(i,lastV)!=Zero) {
                            break;
                        }
                    }
                    IndexType j = min(lastV, prevLastV);
//
//                  T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
//
                    blas::mv(NoTrans, -tau(i),
                             V(_(1,i-1),_(i,j)), V(i,_(i,j)),
                             Zero,
                             _Tr(_(1,i-1),i));
                }
                V(i,i) = Vii;
//
//              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
//
                blas::mv(NoTrans,
                         _Tr(_(1,i-1),_(1,i-1)).upper(),
                         _Tr(_(1,i-1),i));
                Tr(i,i) = tau(i);
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

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename N, typename MV, typename VTAU, typename MT>
void
larft_impl(Direction direction, StoreVectors storeVectors, N n,
           GeMatrix<MV> &V, const DenseVector<VTAU> &tau, TrMatrix<MT> &Tr)
{
    typedef typename TrMatrix<MT>::IndexType  IndexType;

    cxxlapack::larft<IndexType>(getF77Char(direction),
                                getF77Char(storeVectors),
                                n,
                                Tr.dim(),
                                V.data(),
                                V.leadingDimension(),
                                tau.data(),
                                Tr.data(),
                                Tr.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename N, typename MV, typename VTAU, typename MT>
void
larft(Direction direction, StoreVectors storeVectors, N n,
      GeMatrix<MV> &V, const DenseVector<VTAU> &tau, TrMatrix<MT> &Tr)
{
    typedef typename TrMatrix<MT>::IndexType  IndexType;

    LAPACK_DEBUG_OUT("larft");

//
//  Test the input parameters
//
    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(Tr.firstRow()==1);
    ASSERT(Tr.firstCol()==1);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MV>::NoView   V_org = V;
    // copy the full storage!
    typename GeMatrix<MT>::NoView   Tr_org = Tr.general();
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::larft_impl(direction, storeVectors, n, V, tau, Tr);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename GeMatrix<MV>::NoView   V_generic = V;
    typename GeMatrix<MT>::NoView   Tr_generic = Tr.general();

    V = V_org;
    Tr.general() = Tr_org;
//
//  Compare results
//
    external::larft_impl(direction, storeVectors, n, V, tau, Tr);

    bool failed = false;
    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }

    if (! isIdentical(Tr_generic, Tr.general(), "Tr_generic", "_Tr"))
    {
        std::cerr << "CXXLAPACK: Tr_generic = " << Tr_generic << std::endl;
        std::cerr << "F77LAPACK: Tr = " << Tr.general() << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename N, typename MV, typename VTAU, typename MT>
void
larft(Direction direction, StoreVectors storeVectors, N n,
      MV &&V, const VTAU &tau, MT &&Tr)
{
    larft(direction, storeVectors, n, V, tau, Tr);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LARFT_TCC
