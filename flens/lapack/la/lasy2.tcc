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
       SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
      $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LASY2_TCC
#define FLENS_LAPACK_LA_LASY2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename SIGN, typename MTL, typename MTR, typename MB,
          typename SCALE, typename MX, typename XNORM>
typename GeMatrix<MX>::IndexType
lasy2_impl(bool                  transLeft,
           bool                  transRight,
           SIGN                  iSign,
           const GeMatrix<MTL>   &TL,
           const GeMatrix<MTR>   &TR,
           const GeMatrix<MB>    &B,
           SCALE                 &scale,
           GeMatrix<MX>          &X,
           XNORM                 &xNorm)
{
    using std::abs;
    using flens::max;
    using std::swap;

    typedef typename GeMatrix<MX>::ElementType  ElementType;
    typedef typename GeMatrix<MX>::IndexType    IndexType;

    const ElementType Zero(0), Half(0.5), One(1), Two(2), Eight(8);

    const IndexType n1 = TL.numRows();
    const IndexType n2 = TR.numRows();

    const Underscore<IndexType> _;

    IndexType   info = 0;
    IndexType   iPiv;

    bool  bSwap, xSwap;
    ElementType safeMin, beta, gamma, tau1, U11, U12, U22, L21, temp;

//
//    .. Local Arrays ..
//
    IndexType  jPivData_[4],
               locU12Data_[4] = { 3, 4, 1, 2},
               locL21Data_[4] = { 2, 1, 4, 3},
               locU22Data_[4] = { 4, 3, 2, 1};
    DenseVectorView<IndexType>
        jPiv   = typename DenseVectorView<IndexType>::Engine(4, jPivData_),
        locU12 = typename DenseVectorView<IndexType>::Engine(4, locU12Data_),
        locL21 = typename DenseVectorView<IndexType>::Engine(4, locL21Data_),
        locU22 = typename DenseVectorView<IndexType>::Engine(4, locU22Data_);

    ElementType   bTmpData_[4], tmpData_[4], x2Data_[2];
    DenseVectorView<ElementType>
        bTmp = typename DenseVectorView<ElementType>::Engine(4, bTmpData_),
        tmp  = typename DenseVectorView<ElementType>::Engine(4, tmpData_),
        x2   = typename DenseVectorView<ElementType>::Engine(2, x2Data_);

    bool  xSwapPivData_[4] = {false, false, true, true},
          bSwapPivData_[4] = {false, true, false, true};
    DenseVectorView<bool>
        xSwapPiv  = typename DenseVectorView<bool>::Engine(4, xSwapPivData_),
        bSwapPiv  = typename DenseVectorView<bool>::Engine(4, bSwapPivData_);

    ElementType  t16Data_[16];
    GeMatrixView<ElementType>
        T16 = typename GeMatrixView<ElementType>::Engine(4, 4, 4, t16Data_);

//
//  Quick return if possible
//
    if (n1==0 || n2==0) {
        return info;
    }
//
//  Set constants to control overflow
//
    const ElementType eps = lamch<ElementType>(Precision);
    const ElementType smallNum = lamch<ElementType>(SafeMin) / eps;
    const ElementType sign = iSign;

    const IndexType k = n1 + n1 + n2 - 2;

    switch (k) {
//
//  1 by 1: TL11*X + SGN*X*TR11 = B11
//
    case 1:
        tau1 = TL(1,1) + sign*TR(1,1);
        beta = abs(tau1);
        if (beta<=smallNum) {
            tau1 = smallNum;
            beta = smallNum;
            info = 1;
        }

        scale = One;
        gamma = abs(B(1,1));
        if (smallNum*gamma>beta) {
            scale = One/gamma;
        }

        X(1,1) = (B(1,1)*scale) / tau1;
        xNorm = abs(X(1,1));
        return info;

    case 2:
    case 3:
        if (k==2) {
//
//          1 by 2:
//          TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
//                                            [TR21 TR22]
//
            safeMin = max(eps*max(abs(TL(1,1)), abs(TR(1,1)),
                                  abs(TR(1,2)), abs(TR(2,1)),
                                  abs(TR(2,2))),
                          smallNum);
            tmp(1) = TL(1,1) + sign*TR(1,1);
            tmp(4) = TL(1,1) + sign*TR(2,2);
            if (transRight) {
                tmp(2) = sign*TR(2,1);
                tmp(3) = sign*TR(1,2);
            } else {
                tmp(2) = sign*TR(1,2);
                tmp(3) = sign*TR(2,1);
            }
            bTmp(1) = B(1,1);
            bTmp(2) = B(1,2);
        } else {
//
//          2 by 1:
//          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
//            [TL21 TL22] [X21]         [X21]         [B21]
//
            safeMin = max(eps*max(abs(TR(1,1)), abs(TL(1,1)),
                                  abs(TL(1,2)), abs(TL(2,1)),
                                  abs(TL(2,2))),
                          smallNum);
            tmp(1) = TL(1,1) + sign*TR(1,1);
            tmp(4) = TL(2,2) + sign*TR(1,1);
            if (transLeft) {
                tmp(2) = TL(1,2);
                tmp(3) = TL(2,1);
            } else {
                tmp(2) = TL(2,1);
                tmp(3) = TL(1,2);
            }
            bTmp(1) = B(1,1);
            bTmp(2) = B(2,1);
        }
//
//      Solve 2 by 2 system using complete pivoting.
//      Set pivots less than SMIN to SMIN.
//
        iPiv = blas::iamax(tmp);
        U11 = tmp(iPiv);
        if (abs(U11)<=safeMin) {
            info = 1;
            U11 = safeMin;
        }
        U12 = tmp(locU12(iPiv));
        L21 = tmp(locL21(iPiv)) / U11;
        U22 = tmp(locU22(iPiv)) - U12*L21;
        xSwap = xSwapPiv(iPiv);
        bSwap = bSwapPiv(iPiv);
        if (abs(U22)<=safeMin) {
            info = 1;
            U22 = safeMin;
        }
        if (bSwap) {
            temp = bTmp(2);
            bTmp(2) = bTmp(1) - L21*temp;
            bTmp(1) = temp;
        } else {
            bTmp(2) -= L21*bTmp(1);
        }
        scale = One;
        if ((Two*smallNum)*abs(bTmp(2))>abs(U22)
         || (Two*smallNum)*abs(bTmp(1))>abs(U11))
        {
            scale = Half / max(abs(bTmp(1)), abs(bTmp(2)));
            bTmp(1) *= scale;
            bTmp(2) *= scale;
        }
        x2(2) = bTmp(2)/U22;
        x2(1) = bTmp(1)/U11 - (U12/U11)*x2(2);
        if (xSwap) {
            swap(x2(1), x2(2));
        }
        X(1,1) = x2(1);
        if (n1==1) {
            X(1,2) = x2(2);
            xNorm = abs(X(1,1)) + abs(X(1,2));
        } else {
            X(2,1) = x2(2);
            xNorm = max(abs(X(1,1)), abs(X(2,1)));
        }
        return info;

//
//   2 by 2:
//   op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
//     [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
//
//   Solve equivalent 4 by 4 system using complete pivoting.
//   Set pivots less than SMIN to SMIN.
//
    case 4:
        safeMin = max(abs(TR(1,1)), abs(TR(1,2)), abs(TR(2,1)), abs(TR(2,2)));
        safeMin = max(safeMin, abs(TL(1,1)), abs(TL(1,2)),
                               abs(TL(2,1)), abs(TL(2,2)));
        safeMin = max(eps*safeMin, smallNum);
        bTmp(1) = Zero;
        T16 = 0;
        T16(1,1) = TL(1,1) + sign*TR(1,1);
        T16(2,2) = TL(2,2) + sign*TR(1,1);
        T16(3,3) = TL(1,1) + sign*TR(2,2);
        T16(4,4) = TL(2,2) + sign*TR(2,2);
        if (transLeft) {
            T16(1,2) = TL(2,1);
            T16(2,1) = TL(1,2);
            T16(3,4) = TL(2,1);
            T16(4,3) = TL(1,2);
        } else {
            T16(1,2) = TL(1,2);
            T16(2,1) = TL(2,1);
            T16(3,4) = TL(1,2);
            T16(4,3) = TL(2,1);
        }
        if (transRight) {
            T16(1,3) = sign*TR(1,2);
            T16(2,4) = sign*TR(1,2);
            T16(3,1) = sign*TR(2,1);
            T16(4,2) = sign*TR(2,1);
        } else {
            T16(1,3) = sign*TR(2,1);
            T16(2,4) = sign*TR(2,1);
            T16(3,1) = sign*TR(1,2);
            T16(4,2) = sign*TR(1,2);
        }
        bTmp(1) = B(1,1);
        bTmp(2) = B(2,1);
        bTmp(3) = B(1,2);
        bTmp(4) = B(2,2);
//
//      Perform elimination
//
        for (IndexType i=1; i<=3; ++i) {
            ElementType xMax = Zero;
            IndexType ipSv = -1, jpSv = -1;

            for (IndexType ip=i; ip<=4; ++ip) {
                for (IndexType jp=i; jp<=4; ++jp) {
                    if (abs(T16(ip,jp))>=xMax) {
                        xMax = abs(T16(ip,jp));
                        ipSv = ip;
                        jpSv = jp;
                    }
                }
            }
            if (ipSv!=i) {
                blas::swap(T16(ipSv,_), T16(i,_));
                swap(bTmp(i), bTmp(ipSv));
            }
            if (jpSv!=i) {
                blas::swap(T16(_,jpSv), T16(_,i));
            }
            jPiv(i) = jpSv;
            if (abs(T16(i,i))<safeMin) {
                info = 1;
                T16(i,i) = safeMin;
            }
            for (IndexType j=i+1; j<=4; ++j) {
                T16(j,i) /= T16(i,i);
                bTmp(j)  -= T16(j,i)*bTmp(i);
                for (IndexType k=i+1; k<=4; ++k) {
                    T16(j,k) -= T16(j,i)*T16(i,k);
                }
            }
        }
        if (abs(T16(4,4))<safeMin) {
            T16(4,4) = safeMin;
        }
        scale = One;
        if ((Eight*smallNum)*abs(bTmp(1))>abs(T16(1,1))
         || (Eight*smallNum)*abs(bTmp(2))>abs(T16(2,2))
         || (Eight*smallNum)*abs(bTmp(3))>abs(T16(3,3))
         || (Eight*smallNum)*abs(bTmp(4))>abs(T16(4,4)))
        {
            scale = (One/Eight) / max(abs(bTmp(1)), abs(bTmp(2)),
                                      abs(bTmp(3)), abs(bTmp(4)));
            bTmp(1) *= scale;
            bTmp(2) *= scale;
            bTmp(3) *= scale;
            bTmp(4) *= scale;
        }
        for (IndexType i=1; i<=4; ++i) {
            IndexType k = 5 - i;
            const ElementType temp = One/T16(k,k);
            tmp(k) = bTmp(k)*temp;
            for (IndexType j=k+1; j<=4; ++j) {
                tmp(k) -= (temp*T16(k,j))*tmp(j);
            }
        }
        for (IndexType i=1; i<=3; ++i) {
            if (jPiv(4-i)!=4-i) {
                swap(tmp(4-i), tmp(jPiv(4-i)));
            }
        }
        X(1,1) = tmp(1);
        X(2,1) = tmp(2);
        X(1,2) = tmp(3);
        X(2,2) = tmp(4);
        xNorm = max(abs(tmp(1))+abs(tmp(3)), abs(tmp(2))+abs(tmp(4)));
        return info;
    }

    // error if switch does not handle all cases
    ASSERT(0);
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename SIGN, typename MTL, typename MTR, typename MB,
          typename SCALE, typename MX, typename XNORM>
typename GeMatrix<MX>::IndexType
lasy2_impl(bool                  transLeft,
           bool                  transRight,
           SIGN                  sign,
           const GeMatrix<MTL>   &TL,
           const GeMatrix<MTR>   &TR,
           const GeMatrix<MB>    &B,
           SCALE                 &scale,
           GeMatrix<MX>          &X,
           XNORM                 &xNorm)
{
    typedef typename GeMatrix<MX>::IndexType  IndexType;

    IndexType info = cxxlapack::lasy2<IndexType>(transLeft,
                                                 transRight,
                                                 sign,
                                                 TL.numRows(),
                                                 TR.numRows(),
                                                 TL.data(),
                                                 TL.leadingDimension(),
                                                 TR.data(),
                                                 TR.leadingDimension(),
                                                 B.data(),
                                                 B.leadingDimension(),
                                                 scale,
                                                 X.data(),
                                                 X.leadingDimension(),
                                                 xNorm);
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename SIGN, typename MTL, typename MTR, typename MB,
          typename SCALE, typename MX, typename XNORM>
typename GeMatrix<MX>::IndexType
lasy2(bool                  transLeft,
      bool                  transRight,
      SIGN                  sign,
      const GeMatrix<MTL>   &TL,
      const GeMatrix<MTR>   &TR,
      const GeMatrix<MB>    &B,
      SCALE                 &scale,
      GeMatrix<MX>          &X,
      XNORM                 &xNorm)
{
    LAPACK_DEBUG_OUT("lasy2");

    typedef typename GeMatrix<MX>::IndexType IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(sign==1 || sign==-1);

    ASSERT(TL.firstRow()==1);
    ASSERT(TL.firstCol()==1);
    ASSERT(TL.numRows()==TL.numCols());
    ASSERT(TL.numRows()<=2);

    ASSERT(TR.firstRow()==1);
    ASSERT(TR.firstCol()==1);
    ASSERT(TR.numRows()==TR.numCols());
    ASSERT(TR.numRows()<=2);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==TL.numRows());
    ASSERT(B.numCols()==TR.numRows());

    ASSERT(X.firstRow()==1);
    ASSERT(X.firstCol()==1);
    ASSERT(X.numRows()==TL.numRows());
    ASSERT(X.numCols()==TR.numRows());
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    SCALE                           scale_org  = scale;
    typename GeMatrix<MX>::NoView   X_org      = X;
    XNORM                           xNorm_org  = xNorm;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::lasy2_impl(transLeft, transRight, sign,
                                               TL, TR, B,
                                               scale, X, xNorm);
#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    SCALE                           scale_generic   = scale;
    typename GeMatrix<MX>::NoView   X_generic       = X;
    XNORM                           xNorm_generic   = xNorm;

//
//  restore output arguments
//
    scale   = scale_org;
    X       = X_org;
    xNorm   = xNorm_org;

//
//  Compare generic results with results from the native implementation
//

    IndexType info_ = external::lasy2_impl(transLeft, transRight, sign,
                                           TL, TR, B,
                                           scale, X, xNorm);

    bool failed = false;
    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }
    if (! isIdentical(X_generic, X, "X_generic", "X")) {
        std::cerr << "CXXLAPACK: X_generic = "
                  << X_generic << std::endl;
        std::cerr << "F77LAPACK: X = " << X << std::endl;
        failed = true;
    }
    if (! isIdentical(xNorm_generic, xNorm, "xNorm_generic", "xNorm")) {
        std::cerr << "CXXLAPACK: xNorm_generic = "
                  << xNorm_generic << std::endl;
        std::cerr << "F77LAPACK: xNorm = " << xNorm << std::endl;
        failed = true;
    }
    if (! isIdentical(info, info_, " info", "info_")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: info_ = " << info_ << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename SIGN, typename MTL, typename MTR, typename MB,
          typename SCALE, typename MX, typename XNORM>
typename MX::IndexType
lasy2(bool                  transLeft,
      bool                  transRight,
      SIGN                  sign,
      const MTL             &TL,
      const MTR             &TR,
      const MB              &B,
      SCALE                 &&scale,
      MX                    &&X,
      XNORM                 &&xNorm)
{
    typedef typename MX::IndexType IndexType;

    CHECKPOINT_ENTER;
    const IndexType info = lasy2(transLeft, transRight, sign,
                                 TL, TR, B,
                                 scale, X, xNorm);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LASY2_TCC
