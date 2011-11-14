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
      SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
      IMPLICIT NONE
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_AUX_LARFX_TCC
#define FLENS_LAPACK_AUX_LARFX_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larfx_generic(Side side, const DenseVector<VV> &v, const TAU &tau,
              GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::IndexType    IndexType;
    typedef typename GeMatrix<MC>::ElementType  T;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    const T     Zero(0), One(1);

    if (tau==Zero) {
        return;
    }
    if (side==Left) {
//
//      Form  H * C, where H has order m.
//
        switch(m) {
//
//      Special code for 1 x 1 Householder
//
        case 1:
            {
                const T tmp = One - tau*v(1)*v(1);
                for (IndexType j=1; j<=n; ++j) {
                    C(1,j) *= tmp;
                }
            }
            return;
//
//      Special code for 2 x 2 Householder
//
        case 2:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                }
            }
            return;
//
//      Special code for 3 x 3 Householder
//
        case 3:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                }
            }
            return;
//
//      Special code for 4 x 4 Householder
//
        case 4:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                }
            }
            return;
//
//      Special code for 5 x 5 Householder
//
        case 5:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                }
            }
            return;
//
//      Special code for 6 x 6 Householder
//
        case 6:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j) + v6*C(6,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                    C(6,j) -= sum*t6;
                }
            }
            return;
//
//      Special code for 7 x 7 Householder
//
        case 7:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j) + v6*C(6,j) + v7*C(7,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                    C(6,j) -= sum*t6;
                    C(7,j) -= sum*t7;
                }
            }
            return;
//
//      Special code for 8 x 8 Householder
//
        case 8:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                    C(6,j) -= sum*t6;
                    C(7,j) -= sum*t7;
                    C(8,j) -= sum*t8;
                }
            }
            return;
//
//      Special code for 9 x 9 Householder
//
        case 9:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                const T v9 = v(9);
                const T t9 = tau*v9;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j)
                                + v9*C(9,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                    C(6,j) -= sum*t6;
                    C(7,j) -= sum*t7;
                    C(8,j) -= sum*t8;
                    C(9,j) -= sum*t9;
                }
            }
            return;
//
//      Special code for 10 x 10 Householder
//
        case 10:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                const T v9 = v(9);
                const T t9 = tau*v9;
                const T v10 = v(10);
                const T t10 = tau*v10;
                for (IndexType j=1; j<=n; ++j) {
                    const T sum = v1*C(1,j) + v2*C(2,j) + v3*C(3,j) + v4*C(4,j)
                                + v5*C(5,j) + v6*C(6,j) + v7*C(7,j) + v8*C(8,j)
                                + v9*C(9,j) + v10*C(10,j);
                    C(1,j) -= sum*t1;
                    C(2,j) -= sum*t2;
                    C(3,j) -= sum*t3;
                    C(4,j) -= sum*t4;
                    C(5,j) -= sum*t5;
                    C(6,j) -= sum*t6;
                    C(7,j) -= sum*t7;
                    C(8,j) -= sum*t8;
                    C(9,j) -= sum*t9;
                    C(10,j) -= sum*t10;
                }
            }
            return;
//
//      Code for general M
//
        default:
            larf(side, v, tau, C, work);
            return;
        }
    } else {
//
//      Form  C * H, where H has order n.
//
        switch(n) {
//
//      Special code for 1 x 1 Householder
//
        case 1:
            {
                const T tmp = One - tau*v(1)*v(1);
                for (IndexType i=1; i<=m; ++i) {
                    C(i,1) *= tmp;
                }
            }
            return;
//
//      Special code for 2 x 2 Householder
//
        case 2:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                }
            }
            return;
//
//      Special code for 3 x 3 Householder
//
        case 3:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                }
            }
            return;
//
//      Special code for 4 x 4 Householder
//
        case 4:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                }
            }
            return;
//
//      Special code for 5 x 5 Householder
//
        case 5:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                }
            }
            return;
//
//      Special code for 6 x 6 Householder
//
        case 6:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5) + v6*C(i,6);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                    C(i,6) -= sum*t6;
                }
            }
            return;
//
//      Special code for 7 x 7 Householder
//
        case 7:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5) + v6*C(i,6) + v7*C(i,7);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                    C(i,6) -= sum*t6;
                    C(i,7) -= sum*t7;
                }
            }
            return;
//
//      Special code for 8 x 8 Householder
//
        case 8:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5) + v6*C(i,6) + v7*C(i,7) + v8*C(i,8);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                    C(i,6) -= sum*t6;
                    C(i,7) -= sum*t7;
                    C(i,8) -= sum*t8;
                }
            }
            return;
//
//      Special code for 9 x 9 Householder
//
        case 9:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                const T v9 = v(9);
                const T t9 = tau*v9;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5) + v6*C(i,6) + v7*C(i,7) + v8*C(i,8)
                                + v9*C(i,9);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                    C(i,6) -= sum*t6;
                    C(i,7) -= sum*t7;
                    C(i,8) -= sum*t8;
                    C(i,9) -= sum*t9;
                }
            }
            return;
//
//      Special code for 10 x 10 Householder
//
        case 10:
            {
                const T v1 = v(1);
                const T t1 = tau*v1;
                const T v2 = v(2);
                const T t2 = tau*v2;
                const T v3 = v(3);
                const T t3 = tau*v3;
                const T v4 = v(4);
                const T t4 = tau*v4;
                const T v5 = v(5);
                const T t5 = tau*v5;
                const T v6 = v(6);
                const T t6 = tau*v6;
                const T v7 = v(7);
                const T t7 = tau*v7;
                const T v8 = v(8);
                const T t8 = tau*v8;
                const T v9 = v(9);
                const T t9 = tau*v9;
                const T v10 = v(10);
                const T t10 = tau*v10;
                for (IndexType i=1; i<=m; ++i) {
                    const T sum = v1*C(i,1) + v2*C(i,2) + v3*C(i,3) + v4*C(i,4)
                                + v5*C(i,5) + v6*C(i,6) + v7*C(i,7) + v8*C(i,8)
                                + v9*C(i,9) + v10*C(i,10);
                    C(i,1) -= sum*t1;
                    C(i,2) -= sum*t2;
                    C(i,3) -= sum*t3;
                    C(i,4) -= sum*t4;
                    C(i,5) -= sum*t5;
                    C(i,6) -= sum*t6;
                    C(i,7) -= sum*t7;
                    C(i,8) -= sum*t8;
                    C(i,9) -= sum*t9;
                    C(i,10) -= sum*t10;
                }
            }
            return;
//
//      Code for general N
//
        default:
            larf(side, v, tau, C, work);
            return;
        }
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larfx_native(Side side, const DenseVector<VV> &v, const TAU &tau,
             GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::ElementType  T;

    const char      SIDE    = getF77LapackChar<Side>(side);
    const INTEGER   M       = C.numRows();
    const INTEGER   N       = C.numCols();
    const INTEGER   LDC     = C.leadingDimension();

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_IMPL(dlarfx)(&SIDE,
                            &M,
                            &N,
                            v.data(),
                            &tau,
                            C.data(),
                            &LDC,
                            work.data());
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larfx(Side side, const DenseVector<VV> &v, const TAU &tau,
      GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    LAPACK_DEBUG_OUT("larfx");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    ASSERT(v.inc()>0);
    ASSERT(v.firstIndex()==1);
    ASSERT(C.firstRow()==1);
    ASSERT(C.firstCol()==1);

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    if (side==Left) {
        ASSERT(v.length()==m);
        ASSERT(work.length()==n);
    } else {
        ASSERT(v.length()==n);
        ASSERT(work.length()==m);
    }

#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MC>::NoView       C_org      = C;
    typename DenseVector<VV>::NoView    work_org   = work;
#   endif

//
//  Call implementation
//
    larfx_generic(side, v, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MC>::NoView       C_generic      = C;
    typename DenseVector<VV>::NoView    work_generic   = work;

//
//  restore output arguments
//
    C       = C_org;
    work    = work_org;

//
//  Compare generic results with results from the native implementation
//
    larfx_native(side, v, tau, C, work);

    bool failed = false;

    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename VV, typename TAU, typename MC, typename VWORK>
void
larfx(Side side, const VV &v, const TAU &tau, MC &&C, VWORK &&work)
{
    larfx(side, v, tau, C, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFX_TCC
