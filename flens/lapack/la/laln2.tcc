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
       SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
      $                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_LA_LALN2_TCC
#define FLENS_LAPACK_LA_LALN2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename NW, typename SAFEMIN, typename CA, typename MA,
          typename D1, typename D2, typename MB, typename WR, typename WI,
          typename MX, typename SCALE, typename XNORM>
typename GeMatrix<MX>::IndexType
laln2_impl(bool                  transA,
           NW                    nw,
           const SAFEMIN         &safeMin,
           const CA              &ca,
           const GeMatrix<MA>    &A,
           const D1              &d1,
           const D2              &d2,
           const GeMatrix<MB>    &B,
           const WR              &wr,
           const WI              &wi,
           GeMatrix<MX>          &X,
           SCALE                 &scale,
           XNORM                 &xNorm)
{
    using std::abs;
    using std::max;
    using flens::pow;

    typedef typename GeMatrix<MX>::ElementType  T;
    typedef typename GeMatrix<MX>::IndexType    IndexType;

    const T Zero(0), One(1), Two(2);

    const IndexType na = A.numRows();

    IndexType info = 0;

//
//    .. Local Arrays ..
//
    T   _ciData[4], _crData[4];
    GeMatrixView<T> CI = typename GeMatrixView<T>::Engine(2, 2, _ciData, 2);
    GeMatrixView<T> CR = typename GeMatrixView<T>::Engine(2, 2, _crData, 2);

    DenseVectorView<T> civ = typename DenseVectorView<T>::Engine(4, _ciData);
    DenseVectorView<T> crv = typename DenseVectorView<T>::Engine(4, _crData);

    bool _rSwapData[4] = {false, true, false, true};
    bool _zSwapData[4] = {false, false, true, true};
    DenseVectorView<bool>
        rSwap = typename DenseVectorView<bool>::Engine(4, _rSwapData),
        zSwap = typename DenseVectorView<bool>::Engine(4, _zSwapData);

    IndexType _iPivotData[16] = { 1, 2, 3, 4,
                                  2, 1, 4, 3,
                                  3, 4, 1, 2,
                                  4, 3, 2, 1};
    GeMatrixView<IndexType>
        iPivot = typename GeMatrixView<IndexType>::Engine(4, 4, _iPivotData, 4);

//
//  Compute BIGNUM
//
    const T smallNum = Two*lamch<T>(SafeMin);
    const T bigNum = One / smallNum;
    const T safeMini = max(safeMin, smallNum);
//
//  Standard Initializations
//
    scale = One;

    if (na==1) {
//
//      1 x 1  (i.e., scalar) system   C X = B
//
        if (nw==1) {
//
//          Real 1x1 system.
//
//          C = ca A - w D
//
            T csr   = ca*A(1,1) - wr*d1;
            T cNorm = abs(csr);
//
//          If | C | < SMINI, use C = SMINI
//
            if (cNorm<safeMini) {
                csr     = safeMini;
                cNorm   = safeMin;
                info    = 1;
            }
//
//          Check scaling for  X = B / C
//
            T bNorm = abs(B(1,1));
            if (cNorm<One && bNorm>One) {
                if (bNorm>bigNum*cNorm) {
                    scale = One / bNorm;
                }
            }
//
//          Compute X
//
            X(1,1) = (B(1,1)*scale) / csr;
            xNorm = abs(X(1,1));
        } else {
//
//          Complex 1x1 system (w is complex)
//
//          C = ca A - w D
//
            T csr = ca*A(1,1) - wr*d1;
            T csi = -wi*d1;
            T cNorm = abs(csr) + abs(csi);
//
//          If | C | < SMINI, use C = SMINI
//
            if (cNorm<safeMini) {
                csr     = safeMini;
                csi     = Zero;
                cNorm   = safeMini;
                info    = 1;
            }
//
//          Check scaling for  X = B / C
//
            T bNorm = abs(B(1,1)) + abs(B(1,2));
            if (cNorm<One && bNorm>One) {
                if (bNorm>bigNum*cNorm) {
                    scale = One / bNorm;
                }
            }
//
//          Compute X
//
            ladiv(scale*B(1,1), scale*B(1,2), csr, csi, X(1,1), X(1,2));
            xNorm = abs(X(1,1)) + abs(X(1,2));
        }

    } else {
//
//      2x2 System
//
//      Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
//
        CR(1,1) = ca*A(1,1) - wr*d1;
        CR(2,2) = ca*A(2,2) - wr*d2;
        if (transA) {
            CR(1,2) = ca*A(2,1);
            CR(2,1) = ca*A(1,2);
        } else {
            CR(2,1) = ca*A(2,1);
            CR(1,2) = ca*A(1,2);
        }

        if (nw==1) {
//
//          Real 2x2 system  (w is real)
//
//          Find the largest element in C
//
            T cMax = Zero;
            IndexType icMax = 0;

            for (IndexType j=1; j<=4; ++j) {
                if (abs(crv(j))>cMax) {
                    cMax = abs(crv(j));
                    icMax = j;
                }
            }
//
//          If norm(C) < SMINI, use SMINI*identity.
//
            if (cMax<safeMini) {
                const T bNorm = max(abs(B(1,1)), abs(B(2,1)));
                if (safeMini<One && bNorm>One) {
                    if (bNorm>bigNum*safeMini) {
                        scale = One/bNorm;
                    }
                }
                const T temp = scale / safeMini;
                X(1,1) = temp*B(1,1);
                X(2,1) = temp*B(2,1);
                xNorm  = temp*bNorm;
                info   = 1;
                return info;
            }
//
//          Gaussian elimination with complete pivoting.
//
            T UR11  = crv(icMax);
            T CR21  = crv(iPivot(2,icMax));
            T UR12  = crv(iPivot(3,icMax));
            T CR22  = crv(iPivot(4,icMax));
            T UR11R = One / UR11;
            T LR21  = UR11R*CR21;
            T UR22  = CR22 - UR12*LR21;
//
//          If smaller pivot < SMINI, use SMINI
//
            if (abs(UR22)<safeMini) {
                UR22 = safeMini;
                info = 1;
            }
            T BR1, BR2;
            if (rSwap(icMax)) {
                BR1 = B(2,1);
                BR2 = B(1,1);
            } else {
                BR1 = B(1,1);
                BR2 = B(2,1);
            }
            BR2 = BR2 - LR21*BR1;

            const T BBND = max(abs(BR1*(UR22*UR11R)), abs(BR2));
            if (BBND>One && abs(UR22)<One) {
                if (BBND>=bigNum*abs(UR22)) {
                    scale = One / BBND;
                }
            }

            const T XR2 = (BR2*scale) / UR22;
            const T XR1 = (scale*BR1)*UR11R - XR2*(UR11R*UR12);
            if (zSwap(icMax)) {
                X(1,1) = XR2;
                X(2,1) = XR1;
            } else {
                X(1,1) = XR1;
                X(2,1) = XR2;
            }
            xNorm = max(abs(XR1), abs(XR2));
//
//          Further scaling if  norm(A) norm(X) > overflow
//
            if (xNorm>One && cMax>One) {
                if (xNorm>bigNum/cMax) {
                    const T temp = cMax / bigNum;
                    X(1,1) = temp*X(1,1);
                    X(2,1) = temp*X(2,1);
                    xNorm = temp*xNorm;
                    scale = temp*scale;
                }
            }
        } else {
//
//          Complex 2x2 system  (w is complex)
//
//          Find the largest element in C
//
            CI(1,1) = -wi*d1;
            CI(2,1) = Zero;
            CI(1,2) = Zero;
            CI(2,2) = -wi*d2;

            T cMax  = Zero;
            IndexType icMax   = 0;

            for (IndexType j=1; j<=4; ++j) {
                if (abs(crv(j))+abs(civ(j))>cMax) {
                    cMax = abs(crv(j)) + abs(civ(j));
                    icMax = j;
                }
            }
//
//          If norm(C) < SMINI, use SMINI*identity.
//
            if (cMax<safeMini) {
                const T bNorm = max(abs(B(1,1))+abs(B(1,2)),
                                    abs(B(2,1))+abs(B(2,2)));
                if (safeMini<One && bNorm>One) {
                    if (bNorm>bigNum*safeMini) {
                        scale = One / bNorm;
                    }
                }
                const T temp = scale / safeMini;
                X(1,1) = temp*B(1,1);
                X(2,1) = temp*B(2,1);
                X(1,2) = temp*B(1,2);
                X(2,2) = temp*B(2,2);
                xNorm = temp*bNorm;
                info = 1;
                return info;
            }
//
//          Gaussian elimination with complete pivoting.
//
            const T UR11 = crv(icMax);
            const T UI11 = civ(icMax);
            const T CR21 = crv(iPivot(2,icMax));
            const T CI21 = civ(iPivot(2,icMax));
            const T UR12 = crv(iPivot(3,icMax));
            const T UI12 = civ(iPivot(3,icMax));
            const T CR22 = crv(iPivot(4,icMax));
            const T CI22 = civ(iPivot(4,icMax));

            T UR11R, UI11R, LR21, LI21, UR12S, UI12S, UR22, UI22;

            if (icMax==1 || icMax==4) {
//
//              Code when off-diagonals of pivoted C are real
//
                if (abs(UR11)>abs(UI11)) {
                    const T temp = UI11 / UR11;
                    UR11R = One / (UR11*(One+pow(temp,2)));
                    UI11R = -temp*UR11R;
                } else {
                    const T temp = UR11 / UI11;
                    UI11R = -One / (UI11*(One+pow(temp,2)));
                    UR11R = -temp*UI11R;
                }
                LR21 = CR21*UR11R;
                LI21 = CR21*UI11R;
                UR12S = UR12*UR11R;
                UI12S = UR12*UI11R;
                UR22 = CR22 - UR12*LR21;
                UI22 = CI22 - UR12*LI21;
            } else {
//
//              Code when diagonals of pivoted C are real
//
                UR11R = One / UR11;
                UI11R = Zero;
                LR21 = CR21*UR11R;
                LI21 = CI21*UR11R;
                UR12S = UR12*UR11R;
                UI12S = UI12*UR11R;
                UR22 = CR22 - UR12*LR21 + UI12*LI21;
                UI22 = -UR12*LI21 - UI12*LR21;
            }
            const T U22ABS = abs(UR22) + abs(UI22);
//
//          If smaller pivot < SMINI, use SMINI
//
            T BR1, BR2, BI1, BI2;

            if (U22ABS<safeMini) {
                UR22 = safeMini;
                UI22 = Zero;
                info = 1;
            }
            if (rSwap(icMax)) {
                BR2 = B(1,1);
                BR1 = B(2,1);
                BI2 = B(1,2);
                BI1 = B(2,2);
            } else {
                BR1 = B(1,1);
                BR2 = B(2,1);
                BI1 = B(1,2);
                BI2 = B(2,2);
            }
            BR2 = BR2 - LR21*BR1 + LI21*BI1;
            BI2 = BI2 - LI21*BR1 - LR21*BI1;
            const T BBND = max((abs(BR1)+abs(BI1))
                               *(U22ABS*(abs(UR11R)+abs(UI11R))),
                              abs(BR2)+abs(BI2));
            if (BBND>One && U22ABS<One) {
                if (BBND>=bigNum*U22ABS) {
                    scale = One / BBND;
                    BR1 = scale*BR1;
                    BI1 = scale*BI1;
                    BR2 = scale*BR2;
                    BI2 = scale*BI2;
                }
            }

            T XR1, XR2, XI1, XI2;
            ladiv(BR2, BI2, UR22, UI22, XR2, XI2);
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2;
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2;
            if (zSwap(icMax)) {
                X(1,1) = XR2;
                X(2,1) = XR1;
                X(1,2) = XI2;
                X(2,2) = XI1;
            } else {
                X(1,1) = XR1;
                X(2,1) = XR2;
                X(1,2) = XI1;
                X(2,2) = XI2;
            }
            xNorm = max(abs(XR1)+abs(XI1), abs(XR2)+abs(XI2));
//
//          Further scaling if  norm(A) norm(X) > overflow
//
            if (xNorm>One && cMax>One) {
                if (xNorm>bigNum/cMax) {
                    const T temp = cMax / bigNum;
                    X(1,1) = temp*X(1,1);
                    X(2,1) = temp*X(2,1);
                    X(1,2) = temp*X(1,2);
                    X(2,2) = temp*X(2,2);
                    xNorm *= temp;
                    scale *= temp;
                }
            }
        }
    }

    return info;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename NW, typename SAFEMIN, typename CA, typename MA,
          typename D1, typename D2, typename MB, typename WR, typename WI,
          typename MX, typename SCALE, typename XNORM>
typename GeMatrix<MX>::IndexType
laln2_impl(bool                  transA,
           NW                    nw,
           const SAFEMIN         &safeMin,
           const CA              &ca,
           const GeMatrix<MA>    &A,
           const D1              &d1,
           const D2              &d2,
           const GeMatrix<MB>    &B,
           const WR              &wr,
           const WI              &wi,
           GeMatrix<MX>          &X,
           SCALE                 &scale,
           XNORM                 &xNorm)
{
    typedef typename GeMatrix<MX>::IndexType  IndexType;

    IndexType info = cxxlapack::laln2<IndexType>(getF77Char(transA),
                                                 A.numRows(),
                                                 nw,
                                                 safeMin,
                                                 ca,
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 d1,
                                                 d2,
                                                 B.data(),
                                                 B.leadingDimension(),
                                                 wr,
                                                 wi,
                                                 X.data(),
                                                 X.leadingDimension(),
                                                 scale,
                                                 xNorm);
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename NW, typename SAFEMIN, typename CA, typename MA,
          typename D1, typename D2, typename MB, typename WR, typename WI,
          typename MX, typename SCALE, typename XNORM>
typename GeMatrix<MX>::IndexType
laln2(bool                  transA,
      NW                    nw,
      const SAFEMIN         &safeMin,
      const CA              &ca,
      const GeMatrix<MA>    &A,
      const D1              &d1,
      const D2              &d2,
      const GeMatrix<MB>    &B,
      const WR              &wr,
      const WI              &wi,
      GeMatrix<MX>          &X,
      SCALE                 &scale,
      XNORM                 &xNorm)
{
    LAPACK_DEBUG_OUT("BEGIN: laln2");

    typedef typename GeMatrix<MX>::IndexType IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    const IndexType na = A.numRows();
    ASSERT(na==1 || na==2);

    ASSERT(nw==1 || nw==2);

    ASSERT(B.numRows()==na);
    ASSERT(B.numCols()==nw);

    ASSERT(X.numRows()==na);
    ASSERT(X.numCols()==nw);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MX>::NoView   X_org     = X;
    SCALE                           scale_org = scale;
    XNORM                           xNorm_org = xNorm;
#   endif

//
//  Call implementation
//
    const IndexType info = LAPACK_SELECT::laln2_impl(transA, nw, safeMin, ca,
                                                     A, d1, d2, B, wr, wi,
                                                     X, scale, xNorm);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MX>::NoView   X_generic     = X;
    SCALE                           scale_generic = scale;
    XNORM                           xNorm_generic = xNorm;

//
//  restore output arguments
//
    X     = X_org;
    scale = scale_org;
    xNorm = xNorm_org;

//
//  Compare generic results with results from the native implementation
//
    const IndexType _info = external::laln2_impl(transA, nw, safeMin, ca,
                                                 A, d1, d2, B, wr, wi,
                                                 X, scale, xNorm);

    bool failed = false;
    if (! isIdentical(X_generic, X, "X_generic", "X")) {
        std::cerr << "CXXLAPACK: X_generic = " << X_generic << std::endl;
        std::cerr << "F77LAPACK: X = " << X << std::endl;
        failed = true;
    }

    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }

    if (! isIdentical(xNorm_generic, xNorm, "xNorm_generic", "xNorm")) {
        std::cerr << "CXXLAPACK: xNorm_generic = "
                  << xNorm_generic << std::endl;
        std::cerr << "F77LAPACK: xNorm = " << xNorm << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laln2.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: laln2.tcc" << std::endl;
    }
#   endif

    LAPACK_DEBUG_OUT("END: laln2");

    return info;
}


//-- forwarding ----------------------------------------------------------------

template <typename NW, typename SAFEMIN, typename CA, typename MA,
          typename D1, typename D2, typename MB, typename WR, typename WI,
          typename MX, typename SCALE, typename XNORM>
typename MX::IndexType
laln2(bool                  transA,
      NW                    nw,
      const SAFEMIN         &safeMin,
      const CA              &ca,
      const MA              &A,
      const D1              &d1,
      const D2              &d2,
      const MB              &B,
      const WR              &wr,
      const WI              &wi,
      MX                    &&X,
      SCALE                 &scale,
      XNORM                 &xNorm)
{
    return laln2(transA, nw, safeMin, ca, A, d1, d2, B, wr, wi,
                 X, scale, xNorm);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LALN2_TCC
