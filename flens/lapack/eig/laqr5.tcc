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
       SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
      $                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
      $                   LDU, NV, WV, LDWV, NH, WH, LDWH )
 *
 *  -- LAPACK auxiliary routine (version 3.3.0) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2010
 *
 */

#ifndef FLENS_LAPACK_EIG_LAQR5_TCC
#define FLENS_LAPACK_EIG_LAQR5_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename IndexType, typename VSR, typename VSI, typename MH,
          typename MZ, typename MV, typename MU, typename MWV, typename MWH>
void
laqr5_generic(bool                      wantT,
              bool                      wantZ,
              IndexType                 kacc22,
              IndexType                 kTop,
              IndexType                 kBot,
              IndexType                 nShifts,
              DenseVector<VSR>          &sr,
              DenseVector<VSI>          &si,
              GeMatrix<MH>              &H,
              IndexType                 iLoZ,
              IndexType                 iHiZ,
              GeMatrix<MZ>              &Z,
              GeMatrix<MV>              &V,
              GeMatrix<MU>              &U,
              GeMatrix<MWV>             &WV,
              GeMatrix<MWH>             &WH)
{
    using std::abs;
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType  T;

    const T Zero(0), One(1);

    const Underscore<IndexType> _;

    const IndexType n   = H.numRows();
    const IndexType nv  = WV.numRows();
    const IndexType nh  = WH.numCols();

    typedef typename GeMatrix<MH>::VectorView   VectorView;
    T           vtBuffer[3];
    VectorView  vt = typename VectorView::Engine(3, vtBuffer);

//
//  ==== If there are no shifts, then there is nothing to do. ====
//
    if (nShifts<2) {
        return;
    }
//
//  ==== If the active block is empty or 1-by-1, then there
//  .    is nothing to do. ====
//
    if (kTop>=kBot) {
        return;
    }
//
//  ==== Shuffle shifts into pairs of real shifts and pairs
//  .    of complex conjugate shifts assuming complex
//  .    conjugate shifts are already adjacent to one
//  .    another. ====
//
    for (IndexType i=1; i<=nShifts-2; i+=2) {
        if (si(i)!=-si(i+1)) {
            T tmp = sr(i);
            sr(i) = sr(i+1);
            sr(i+1) = sr(i+2);
            sr(i+2) = tmp;

            tmp = si(i);
            si(i) = si(i+1);
            si(i+1) = si(i+2);
            si(i+2) = tmp;
        }
    }
//
//  ==== NSHFTS is supposed to be even, but if it is odd,
//  .    then simply reduce it by one.  The shuffle above
//  .    ensures that the dropped shift is real and that
//  .    the remaining shifts are paired. ====
//
    const IndexType ns = nShifts - (nShifts % 2);
//
//  ==== Machine constants for deflation ====
//
    T safeMin   = lamch<T>(SafeMin);
    T safeMax   = One/safeMin;
    labad(safeMin, safeMax);

    const T ulp         = lamch<T>(Precision);
    const T smallNum    = safeMin*(T(n)/ulp);
//
//  ==== Use accumulated reflections to update far-from-diagonal
//  .    entries ? ====
//
    const bool accum = (kacc22==1) || (kacc22==2);
//
//  ==== If so, exploit the 2-by-2 block structure? ====
//
    const bool blk22 = (ns>2) && (kacc22==2);
//
//  ==== clear trash ====
//
    if (kTop+2<=kBot) {
        H(kTop+2,kTop) = Zero;
    }
//
//  ==== NBMPS = number of 2-shift bulges in the chain ====
//
    const IndexType nBmps = ns/2;
//
//  ==== KDU = width of slab ====
//
    const IndexType kdu = 6*nBmps - 3;
//
//  ==== Create and chase chains of NBMPS bulges ====
//
    for (IndexType inCol=3*(1-nBmps)+kTop-1; inCol<=kBot-2; inCol+=3*nBmps-2) {
        IndexType ndCol = inCol + kdu;
        if (accum) {
            auto U_ = U(_(1,kdu),_(1,kdu));
            U_           = Zero;
            U_.diag(0)   = One;
        }
//
//      ==== Near-the-diagonal bulge chase.  The following loop
//      .    performs the near-the-diagonal part of a small bulge
//      .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
//      .    chunk extends from column INCOL to column NDCOL
//      .    (including both column INCOL and column NDCOL). The
//      .    following loop chases a 3*NBMPS column long chain of
//      .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
//      .    may be less than KTOP and and NDCOL may be greater than
//      .    KBOT indicating phantom columns from which to chase
//      .    bulges before they are actually introduced or to which
//      .    to chase bulges beyond column KBOT.)  ====
//
        const IndexType krColMax = min(inCol+3*nBmps-3, kBot-2);

        for (IndexType krCol=inCol; krCol<=krColMax; ++krCol) {
//
//          ==== Bulges number MTOP to MBOT are active double implicit
//          .    shift bulges.  There may or may not also be small
//          .    2-by-2 bulge, if there is room.  The inactive bulges
//          .    (if any) must wait until the active bulges have moved
//          .    down the diagonal to make room.  The phantom matrix
//          .    paradigm described above helps keep track.  ====
//
            const IndexType mTop = max(IndexType(1), ((kTop-1)-krCol+2)/3+1);
            const IndexType mBot = min(nBmps, (kBot-krCol)/3);
            const IndexType m22  = mBot + 1;
            const bool bmp22     = (mBot<nBmps ) && (krCol+3*(m22-1))==(kBot-2);
//
//          ==== Generate reflections to chase the chain right
//          .    one column.  (The minimum value of K is KTOP-1.) ====
//
            IndexType   k;
            T           alpha, beta;

            for (IndexType m=mTop; m<=mBot; ++m) {
                k = krCol + 3*(m-1);
                if (k==kTop-1) {
                    laqr1(H(_(kTop,kTop+2),_(kTop,kTop+2)),
                          sr(2*m-1), si(2*m-1), sr(2*m), si(2*m),
                          V(_(1,3),m));
                    alpha = V(1,m);
                    larfg(IndexType(3), alpha, V(_(2,3),m), V(1,m));
                } else {
                    beta = H(k+1,k);
                    V(2,m) = H(k+2,k);
                    V(3,m) = H(k+3,k);
                    larfg(IndexType(3), beta, V(_(2,3),m), V(1,m));
//
//                  ==== A Bulge may collapse because of vigilant
//                  .    deflation or destructive underflow.  In the
//                  .    underflow case, try the two-small-subdiagonals
//                  .    trick to try to reinflate the bulge.  ====
//
                    if (H(k+3,k)!=Zero
                     || H(k+3,k+1)!=Zero
                     || H(k+3,k+2)==Zero) {
//
//                   ==== Typical case: not collapsed (yet). ====
//
                        H(k+1, k) = beta;
                        H(k+2, k) = Zero;
                        H(k+3, k) = Zero;
                    } else {
//
//                      ==== Atypical case: collapsed.  Attempt to
//                      .    reintroduce ignoring H(K+1,K) and H(K+2,K).
//                      .    If the fill resulting from the new
//                      .    reflector is too large, then abandon it.
//                      .    Otherwise, use the new one. ====
//
                        laqr1(H(_(k+1,k+3),_(k+1,k+3)),
                              sr(2*m-1), si(2*m-1), sr(2*m), si(2*m),
                              vt);
                        alpha = vt(1);
                        larfg(3, alpha, vt(_(2,3)), vt(1));
                        const T refSum = vt(1)*(H(k+1,k) + vt(2)*H(k+2,k));
                        if (abs(H(k+2,k)-refSum*vt(2)) + abs(refSum*vt(3))
                            > ulp*(abs(H(k,k))+abs(H(k+1,k+1))+abs(H(k+2,k+2))))
                        {
//
//                          ==== Starting a new bulge here would
//                          .    create non-negligible fill.  Use
//                          .    the old one with trepidation. ====
//
                            H(k+1, k) = beta;
                            H(k+2, k) = Zero;
                            H(k+3, k) = Zero;
                        } else {
//
//                          ==== Stating a new bulge here would
//                          .    create only negligible fill.
//                          .    Replace the old reflector with
//                          .    the new one. ====
//
                            H(k+1, k) -= refSum;
                            H(k+2, k) = Zero;
                            H(k+3, k) = Zero;
                            V(1, m) = vt(1);
                            V(2, m) = vt(2);
                            V(3, m) = vt(3);
                        }
                    }
                }
            }
//
//          ==== Generate a 2-by-2 reflection, if needed. ====
//
            k = krCol + 3*(m22-1);
            if (bmp22) {
                if (k==kTop-1) {
                    laqr1(H(_(k+1,k+2),_(k+1,k+2)),
                          sr(2*m22-1), si(2*m22-1), sr(2*m22), si(2*m22),
                          V(_(1,2),m22));
                    beta = V(1, m22);
                    larfg(IndexType(2), beta, V(_(2,2),m22), V(1,m22));
                } else {
                    beta = H(k+1, k);
                    V(2, m22) = H(k+2,k);
                    larfg(IndexType(2), beta, V(_(2,2),m22), V(1,m22));
                    H(k+1, k) = beta;
                    H(k+2, k) = Zero;
                }
            }
//
//          ==== Multiply H by reflections from the left ====
//
            IndexType jBot;
            if (accum) {
                jBot = min(ndCol, kBot);
            } else if (wantT) {
                jBot = n;
            } else {
                jBot = kBot;
            }
            for (IndexType j=max(kTop,krCol); j<=jBot; ++j) {
                IndexType mEnd = min(mBot, (j-krCol+2)/3);
                for (IndexType m=mTop; m<=mEnd; ++m) {
                    k = krCol + 3*(m-1);
                    const T refSum = V(1,m)*(H(k+1,j) + V(2,m)*H(k+2,j)
                                                      + V(3,m)*H(k+3,j));
                    H(k+1,j) -= refSum;
                    H(k+2,j) -= refSum*V(2,m);
                    H(k+3,j) -= refSum*V(3,m);
                }
            }
            if (bmp22) {
                k = krCol + 3*(m22-1);
                for (IndexType j=max(k+1,kTop); j<=jBot; ++j) {
                    const T refSum = V(1,m22)*(H(k+1,j)+V(2,m22)*H(k+2,j));
                    H(k+1,j) -= refSum;
                    H(k+2,j) -= refSum*V(2,m22);
                }
            }
//
//          ==== Multiply H by reflections from the right.
//          .    Delay filling in the last row until the
//          .    vigilant deflation check is complete. ====
//
            IndexType jTop;
            if (accum) {
                jTop = max(kTop, inCol);
            } else if (wantT) {
                jTop = 1;
            } else {
                jTop = kTop;
            }
            for (IndexType m=mTop; m<=mBot; ++m) {
                if (V(1,m)!=Zero) {
                    k = krCol + 3*(m-1);
                    for (IndexType j=jTop; j<=min(kBot,k+3); ++j) {
                        const T refSum = V(1,m)*(H(j,k+1) + V(2,m)*H(j,k+2)
                                                          + V(3,m)*H(j,k+3));
                        H(j, k+1) -= refSum;
                        H(j, k+2) -= refSum*V(2, m);
                        H(j, k+3) -= refSum*V(3, m);
                    }

                    if (accum) {
//
//                   ==== Accumulate U. (If necessary, update Z later
//                   .    with with an efficient matrix-matrix
//                   .    multiply.) ====
//
                        IndexType kms   = k - inCol;
                        IndexType j1    = max(IndexType(1),kTop-inCol);

                        for (IndexType j=j1; j<=kdu; ++j) {
                            const T refSum = V(1,m)*(U(j,kms+1)
                                                        + V(2,m)*U(j,kms+2)
                                                        + V(3,m)*U(j,kms+3));
                            U(j, kms+1) -= refSum;
                            U(j, kms+2) -= refSum*V(2, m);
                            U(j, kms+3) -= refSum*V(3, m);
                        }
                    } else if (wantZ) {
//
//                      ==== U is not accumulated, so update Z
//                      .    now by multiplying by reflections
//                      .    from the right. ====
//
                        for (IndexType j=iLoZ; j<=iHiZ; ++j) {
                            const T refSum = V(1,m)*(Z(j,k+1)
                                                        +V(2,m)*Z(j,k+2)
                                                        +V(3,m)*Z(j,k+3));
                            Z(j,k+1) -= refSum;
                            Z(j,k+2) -= refSum*V(2,m);
                            Z(j,k+3) -= refSum*V(3,m);
                        }
                    }
                }
            }
//
//          ==== Special case: 2-by-2 reflection (if needed) ====
//
            k = krCol + 3*(m22-1);
            if (bmp22) {
                if (V(1,m22)!=Zero) {
                    for (IndexType j=jTop; j<=min(kBot,k+3); ++j) {
                        const T refSum = V(1,m22)*(H(j,k+1)+V(2,m22)*H(j,k+2));
                        H(j,k+1) -= refSum;
                        H(j,k+2) -= refSum*V(2,m22);
                    }

                    if (accum) {
                        IndexType kms   = k - inCol;
                        IndexType j1    = max(IndexType(1),kTop-inCol);

                        for (IndexType j=j1; j<=kdu; ++j) {
                            const T refSum = V(1,m22)*(U(j,kms+1)
                                                        + V(2,m22)*U(j,kms+2));
                            U(j,kms+1 ) -= refSum;
                            U(j,kms+2 ) -= refSum*V(2,m22);
                        }
                    } else if (wantZ) {
                        for (IndexType j=iLoZ; j<=iHiZ; ++j) {
                            const T refSum = V(1,m22)*(Z(j,k+1)
                                                        + V(2,m22)*Z(j,k+2));
                            Z(j,k+1) -= refSum;
                            Z(j,k+2) -= refSum*V(2,m22);
                        }
                    }
                }
            }
//
//          ==== Vigilant deflation check ====
//
            IndexType mStart = mTop;
            if (krCol+3*(mStart-1)<kTop) {
                ++mStart;
            }
            IndexType mEnd = mBot;
            if (bmp22) {
                ++mEnd;
            }
            if (krCol==kBot-2) {
                ++mEnd;
            }
            for (IndexType m=mStart; m<=mEnd; ++m) {
                k = min(kBot-1, krCol+3*(m-1));
//
//              ==== The following convergence test requires that
//              .    the tradition small-compared-to-nearby-diagonals
//              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
//              .    criteria both be satisfied.  The latter improves
//              .    accuracy in some examples. Falling back on an
//              .    alternate convergence criterion when TST1 or TST2
//              .    is zero (as done here) is traditional but probably
//              .    unnecessary. ====
//
                if (H(k+1,k)!=Zero) {
                    T test1 = abs(H(k,k)) + abs(H(k+1,k+1));
                    if (test1==Zero) {
                        if (k>=kTop+1) {
                            test1 += abs(H(k,k-1));
                        }
                        if (k>=kTop+2) {
                            test1 += abs(H(k,k-2));
                        }
                        if (k>=kTop+3) {
                            test1 += abs(H(k,k-3));
                        }
                        if (k<=kBot-2) {
                            test1 += abs(H(k+2,k+1));
                        }
                        if (k<=kBot-3) {
                            test1 += abs(H(k+3,k+1));
                        }
                        if (k<=kBot-4) {
                            test1 += abs(H(k+4,k+1));
                        }
                    }
                    if (abs(H(k+1,k))<=max(smallNum, ulp*test1)) {
                        const T H12 = max(abs(H(k+1,k)), abs(H(k,k+1)));
                        const T H21 = min(abs(H(k+1,k)), abs(H(k,k+1)));
                        const T H11 = max(abs(H(k+1,k+1)),
                                          abs(H(k,k)-H(k+1,k+1)));
                        const T H22 = min(abs(H(k+1,k+1)),
                                          abs(H(k,k)-H(k+1,k+1)));
                        const T scal = H11 + H12;
                        const T test2 = H22*(H11/scal);

                        if (test2==Zero
                         || H21*(H12/scal)<=max(smallNum,ulp*test2))
                        {
                            H(k+1,k) = Zero;
                        }
                    }
                }
            }
//
//          ==== Fill in the last row of each bulge. ====
//
            mEnd = min(nBmps, (kBot-krCol-1)/3);
            for (IndexType m=mTop; m<=mEnd; ++m) {
                k = krCol + 3*(m-1);
                const T refSum = V(1,m)*V(3,m)*H(k+4,k+3);
                H(k+4,k+1) = -refSum;
                H(k+4,k+2) = -refSum*V(2,m);
                H(k+4,k+3) -= refSum*V(3,m);
            }
//
//          ==== End of near-the-diagonal bulge chase. ====
//
        }
//
//      ==== Use U (if accumulated) to update far-from-diagonal
//      .    entries in H.  If required, use U to update Z as
//      .    well. ====
//
        if (accum) {
            IndexType jTop, jBot;
            if (wantT) {
                jTop = 1;
                jBot = n;
            } else {
                jTop = kTop;
                jBot = kBot;
            }
            if ((!blk22) || (inCol<kTop) || (ndCol>kBot) || (ns<=2)) {
//
//              ==== Updates not exploiting the 2-by-2 block
//              .    structure of U.  K1 and NU keep track of
//              .    the location and size of U in the special
//              .    cases of introducing bulges and chasing
//              .    bulges off the bottom.  In these special
//              .    cases and in case the number of shifts
//              .    is NS = 2, there is no 2-by-2 block
//              .    structure to exploit.  ====
//
                const IndexType k1 = max(IndexType(1), kTop-inCol);
                const IndexType nu  = (kdu-max(IndexType(0), ndCol-kBot)) -k1+1;
                const IndexType _nu = (kdu-max(IndexType(0), ndCol-kBot));
//
//              ==== Horizontal Multiply ====
//
                for (IndexType jCol=min(ndCol,kBot)+1; jCol<=jBot; jCol+=nh) {
                    const IndexType jLen = min(nh, jBot-jCol+1);

                    auto _U     = U(_(k1, _nu), _(k1, _nu));
                    auto _H     = H(_(inCol+k1,inCol+_nu),_(jCol,jCol+jLen-1));
                    auto _WH    = WH(_(1,nu),_(1,jLen));

                    blas::mm(ConjTrans, NoTrans, One, _U, _H, Zero, _WH);
                    _H = _WH;
                }
//
//              ==== Vertical multiply ====
//
                for (IndexType jRow=jTop; jRow<=max(kTop,inCol)-1; jRow+=nv) {
                    const IndexType jLen = min(nv, max(kTop,inCol)-jRow);

                    auto _H     = H(_(jRow,jRow+jLen-1),_(inCol+k1,inCol+_nu));
                    auto _U     = U(_(k1,_nu),_(k1,_nu));
                    auto _WV    = WV(_(1,jLen),_(1,nu));

                    blas::mm(NoTrans, NoTrans, One, _H, _U, Zero, _WV);
                    _H = _WV;
                }
//
//              ==== Z multiply (also vertical) ====
//
                if (wantZ) {
                    for (IndexType jRow=iLoZ; jRow<=iHiZ; jRow+=nv) {
                        const IndexType jLen = min(nv, iHiZ-jRow+1);

                        auto _Z  = Z(_(jRow,jRow+jLen-1),_(inCol+k1,inCol+_nu));
                        auto _U  = U(_(k1,_nu),_(k1,_nu));
                        auto _WV = WV(_(1,jLen),_(1,nu));

                        blas::mm(NoTrans, NoTrans, One, _Z, _U, Zero, _WV);
                        _Z = _WV;
                    }
                }
            } else {
//
//              ==== Updates exploiting U's 2-by-2 block structure.
//              .    (I2, I4, J2, J4 are the last rows and columns
//              .    of the blocks.) ====
//
                const IndexType i2 = (kdu+1 ) / 2;
                const IndexType i4 = kdu;
                const IndexType j2 = i4 - i2;
                const IndexType j4 = kdu;
//
//              ==== KZS and KNZ deal with the band of zeros
//              .    along the diagonal of one of the triangular
//              .    blocks. ====
//
                const IndexType kZs = (j4-j2) - (ns+1);
                const IndexType kNz = ns + 1;
//
//              ==== Horizontal multiply ====
//
                for (IndexType jCol=min(ndCol, kBot)+1; jCol<=jBot; jCol+=nh) {
                    const IndexType jLen = min(nh, jBot-jCol+1);
//
//                  ==== Copy bottom of H to top+KZS of scratch ====
//                  (The first KZS rows get multiplied by zero.) ====
//
                    WH(_(kZs+1,kZs+kNz),_(1,jLen))
                        = H(_(inCol+j2+1,inCol+j2+kNz),_(jCol,jCol+jLen-1));
//
//                  ==== Multiply by U21**T ====
//
                    WH(_(1,kZs),_(1,jLen)) = Zero;
                    blas::mm(Left, ConjTrans, One,
                             U(_(j2+1,j2+kNz),_(kZs+1,kZs+kNz)).upper(),
                             WH(_(kZs+1,kZs+kNz),_(1,jLen)));
//
//                  ==== Multiply top of H by U11**T ====
//
                    blas::mm(ConjTrans, NoTrans, One,
                             U(_(1,j2),_(1,i2)),
                             H(_(inCol+1,inCol+j2),_(jCol,jCol+jLen-1)),
                             One,
                             WH(_(1,i2),_(1,jLen)));
//
//                  ==== Copy top of H to bottom of WH ====
//
                    WH(_(i2+1,i2+j2),_(1,jLen))
                        = H(_(inCol+1,inCol+j2),_(jCol,jCol+jLen-1));
//
//                  ==== Multiply by U21**T ====
//
                    blas::mm(Left, ConjTrans, One,
                             U(_(1,j2),_(i2+1,i2+j2)).lower(),
                             WH(_(i2+1,i2+j2),_(1,jLen)));
//
//                  ==== Multiply by U22 ====
//
                    blas::mm(ConjTrans, NoTrans, One,
                             U(_(j2+1,j4),_(i2+1,i4)),
                             H(_(inCol+j2+1,inCol+j4),_(jCol,jCol+jLen-1)),
                             One,
                             WH(_(i2+1,i4),_(1,jLen)));
//
//                  ==== Copy it back ====
//
                    H(_(inCol+1,inCol+kdu),_(jCol,jCol+jLen-1))
                        = WH(_(1,kdu),_(1,jLen));
                }
//
//              ==== Vertical multiply ====
//
                for (IndexType jRow=jTop; jRow<=max(inCol,kTop)-1; jRow+=nv) {
                    const IndexType jLen = min(nv, max(inCol,kTop) -jRow);
//
//                  ==== Copy right of H to scratch (the first KZS
//                  .    columns get multiplied by zero) ====
//
                    WV(_(1,jLen),_(kZs+1,kZs+kNz))
                        = H(_(jRow,jRow+jLen-1),_(inCol+j2+1, inCol+j2+kNz));
//
//                  ==== Multiply by U21 ====
//
                    WV(_(1,jLen),_(1,kZs)) = Zero;
                    blas::mm(Right, NoTrans, One,
                             U(_(j2+1,j2+kNz),_(kZs+1,kZs+kNz)).upper(),
                             WV(_(1,jLen),_(kZs+1,kZs+kNz)));
//
//                  ==== Multiply by U11 ====
//
                    blas::mm(NoTrans, NoTrans, One,
                             H(_(jRow,jRow+jLen-1),_(inCol+1,inCol+j2)),
                             U(_(1,j2),_(1,i2)),
                             One,
                             WV(_(1,jLen),_(1,i2)));
//
//                  ==== Copy left of H to right of scratch ====
//
                    WV(_(1,jLen),_(i2+1,i2+j2))
                        = H(_(jRow,jRow+jLen-1), _(inCol+1,inCol+j2));
//
//                  ==== Multiply by U21 ====
//
                    blas::mm(Right, NoTrans, One,
                             U(_(1,i4-i2),_(i2+1,i4)).lower(),
                             WV(_(1,jLen),_(i2+1,i4)));
//
//                  ==== Multiply by U22 ====
//
                    blas::mm(NoTrans, NoTrans, One,
                             H(_(jRow,jRow+jLen-1),_(inCol+j2+1,inCol+j4)),
                             U(_(j2+1,j4),_(i2+1,i4)),
                             One,
                             WV(_(1,jLen),_(i2+1,i4)));
//
//                  ==== Copy it back ====
//
                    H(_(jRow,jRow+jLen-1),_(inCol+1,inCol+kdu))
                        = WV(_(1,jLen),_(1,kdu));
                }
//
//              ==== Multiply Z (also vertical) ====
//
                if (wantZ) {
                    for (IndexType jRow=iLoZ; jRow<=iHiZ; jRow+=nv) {
                        const IndexType jLen = min(nv,iHiZ-jRow+1);
//
//                      ==== Copy right of Z to left of scratch (first
//                      .     KZS columns get multiplied by zero) ====
//
                        WV(_(1,jLen),_(kZs+1,kZs+kNz))
                            = Z(_(jRow,jRow+jLen-1),_(inCol+j2+1,inCol+j2+kNz));
//
//                      ==== Multiply by U12 ====
//
                        WV(_(1,jLen),_(1,kZs)) = Zero;
                        blas::mm(Right, NoTrans, One,
                                 U(_(j2+1,j2+kNz), _(kZs+1,kZs+kNz)).upper(),
                                 WV(_(1,jLen),_(kZs+1,kZs+kNz)));
//
//                      ==== Multiply by U11 ====
//
                        blas::mm(NoTrans, NoTrans, One,
                                 Z(_(jRow,jRow+jLen-1), _(inCol+1,inCol+j2)),
                                 U(_(1,j2),_(1,i2)),
                                 One,
                                 WV(_(1,jLen),_(1,i2)));
//
//                      ==== Copy left of Z to right of scratch ====
//
                        WV(_(1,jLen),_(i2+1,i2+j2))
                            = Z(_(jRow,jRow+jLen-1), _(inCol+1,inCol+j2));
//
//                      ==== Multiply by U21 ====
//
                        blas::mm(Right, NoTrans, One,
                                 U(_(1,i4-i2),_(i2+1,i4)).lower(),
                                 WV(_(1,jLen),_(i2+1,i4)));
//
//                      ==== Multiply by U22 ====
//
                        blas::mm(NoTrans, NoTrans, One,
                                 Z(_(jRow,jRow+jLen-1),_(inCol+j2+1,inCol+j4)),
                                 U(_(j2+1,j4),_(i2+1,i4)),
                                 One,
                                 WV(_(1,jLen),_(i2+1,i4)));
//
//                      ==== Copy the result back to Z ====
//
                        Z(_(jRow,jRow+jLen-1),_(inCol+1,inCol+kdu))
                            = WV(_(1,jLen),_(1,kdu));
                    }
                }
            }
        }
    }
}

//== interface for native lapack ===============================================

#ifdef TODO_CHECK_CXXLAPACK

template <typename IndexType, typename VSR, typename VSI, typename MH,
          typename MZ, typename MV, typename MU, typename MWV, typename MWH>
void
laqr5_native(bool                      wantT,
             bool                      wantZ,
             IndexType                 kacc22,
             IndexType                 kTop,
             IndexType                 kBot,
             IndexType                 nShifts,
             DenseVector<VSR>          &sr,
             DenseVector<VSI>          &si,
             GeMatrix<MH>              &H,
             IndexType                 iLoZ,
             IndexType                 iHiZ,
             GeMatrix<MZ>              &Z,
             GeMatrix<MV>              &V,
             GeMatrix<MU>              &U,
             GeMatrix<MWV>             &WV,
             GeMatrix<MWH>             &WH)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    const LOGICAL    WANTT      = wantT;
    const LOGICAL    WANTZ      = wantZ;
    const INTEGER    KACC22     = kacc22;
    const INTEGER    N          = H.numRows();
    const INTEGER    KTOP       = kTop;
    const INTEGER    KBOT       = kBot;
    const INTEGER    NSHFTS     = nShifts;
    const INTEGER    LDH        = H.leadingDimension();
    const INTEGER    ILOZ       = iLoZ;
    const INTEGER    IHIZ       = iHiZ;
    const INTEGER    LDZ        = Z.leadingDimension();
    const INTEGER    LDV        = V.leadingDimension();
    const INTEGER    LDU        = U.leadingDimension();
    const INTEGER    NV         = WV.numRows();
    const INTEGER    LDWV       = WV.leadingDimension();
    const INTEGER    NH         = WH.numCols();
    const INTEGER    LDWH       = WH.leadingDimension();

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_IMPL(dlaqr5)(&WANTT,
                            &WANTZ,
                            &KACC22,
                            &N,
                            &KTOP,
                            &KBOT,
                            &NSHFTS,
                            sr.data(),
                            si.data(),
                            H.data(),
                            &LDH,
                            &ILOZ,
                            &IHIZ,
                            Z.data(),
                            &LDZ,
                            V.data(),
                            &LDV,
                            U.data(),
                            &LDU,
                            &NV,
                            WV.data(),
                            &LDWV,
                            &NH,
                            WH.data(),
                            &LDWH);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================
template <typename IndexType, typename VSR, typename VSI, typename MH,
          typename MZ, typename MV, typename MU, typename MWV, typename MWH>
void
laqr5(bool                      wantT,
      bool                      wantZ,
      IndexType                 kacc22,
      IndexType                 kTop,
      IndexType                 kBot,
      IndexType                 nShifts,
      DenseVector<VSR>          &sr,
      DenseVector<VSI>          &si,
      GeMatrix<MH>              &H,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      GeMatrix<MZ>              &Z,
      GeMatrix<MV>              &V,
      GeMatrix<MU>              &U,
      GeMatrix<MWV>             &WV,
      GeMatrix<MWH>             &WH)
{
    LAPACK_DEBUG_OUT("laqr5");

    using std::max;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT((kacc22==0)||(kacc22==1)||(kacc22==2));
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    ASSERT(1<=kTop);
    ASSERT(kBot<=n);

    ASSERT(nShifts>0);
    ASSERT(nShifts % 2 == 0);

    ASSERT(sr.length()==nShifts);
    ASSERT(si.length()==nShifts);

    if (wantZ) {
        ASSERT(1<=iLoZ);
        ASSERT(iLoZ<=iHiZ);
        ASSERT(iHiZ<=n);
    }

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);
    ASSERT(V.numRows()>=3);
    ASSERT(V.numCols()==nShifts/2);

    ASSERT(U.firstRow()==1);
    ASSERT(U.firstCol()==1);
    ASSERT(U.numRows()>=3*nShifts-3);

    ASSERT(WH.firstRow()==1);
    ASSERT(WH.firstCol()==1);
    ASSERT(WH.numRows()>=3*nShifts-3);
    ASSERT(WH.numCols()>=1);

    ASSERT(WV.firstRow()==1);
    ASSERT(WV.firstCol()==1);
    ASSERT(WV.numRows()>=1);
    ASSERT(WV.numCols()>=3*nShifts-3);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename DenseVector<VSR>::NoView   sr_org = sr;
    typename DenseVector<VSI>::NoView   si_org = si;
    typename GeMatrix<MH>::NoView       H_org  = H;
    typename GeMatrix<MZ>::NoView       Z_org  = Z;
    typename GeMatrix<MV>::NoView       V_org  = V;
    typename GeMatrix<MU>::NoView       U_org  = U;
    typename GeMatrix<MWV>::NoView      WV_org = WV;
    typename GeMatrix<MWH>::NoView      WH_org = WH;
#   endif

//
//  Call implementation
//
    laqr5_generic(wantT, wantZ, kacc22, kTop, kBot, nShifts,
                  sr, si, H, iLoZ, iHiZ, Z,
                  V, U, WV, WH);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename DenseVector<VSR>::NoView   sr_generic = sr;
    typename DenseVector<VSI>::NoView   si_generic = si;
    typename GeMatrix<MH>::NoView       H_generic  = H;
    typename GeMatrix<MZ>::NoView       Z_generic  = Z;
    typename GeMatrix<MV>::NoView       V_generic  = V;
    typename GeMatrix<MU>::NoView       U_generic  = U;
    typename GeMatrix<MWV>::NoView      WV_generic = WV;
    typename GeMatrix<MWH>::NoView      WH_generic = WH;
//
//  restore output arguments
//
    sr = sr_org;
    si = si_org;
    H  = H_org;
    Z  = Z_org;
    V  = V_org;
    U  = U_org;
    WV = WV_org;
    WH = WH_org;

//
//  Compare generic results with results from the native implementation
//
    laqr5_native(wantT, wantZ, kacc22, kTop, kBot, nShifts,
                 sr, si, H, iLoZ, iHiZ, Z,
                 V, U, WV, WH);

    bool failed = false;
    if (! isIdentical(sr_generic, sr, "sr_generic", "sr")) {
//        std::cerr << "CXXLAPACK: sr_generic = " << sr_generic << std::endl;
//        std::cerr << "F77LAPACK: sr = " << sr << std::endl;
        failed = true;
    }

    if (! isIdentical(si_generic, si, "si_generic", "si")) {
//        std::cerr << "CXXLAPACK: si_generic = " << si_generic << std::endl;
//        std::cerr << "F77LAPACK: si = " << si << std::endl;
        failed = true;
    }

    if (! isIdentical(H_generic, H, "H_generic", "H")) {
//        std::cerr << "CXXLAPACK: H_generic = " << H_generic << std::endl;
//        std::cerr << "F77LAPACK: H = " << H << std::endl;
        failed = true;
    }

    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
//        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
//        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
        failed = true;
    }

    if (! isIdentical(V_generic, V, "V_generic", "V")) {
//        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
//        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }

    if (! isIdentical(U_generic, U, "U_generic", "U")) {
//        std::cerr << "CXXLAPACK: U_generic = " << U_generic << std::endl;
//        std::cerr << "F77LAPACK: U = " << U << std::endl;
        failed = true;
    }

    if (! isIdentical(WV_generic, WV, "WV_generic", "WV")) {
//        std::cerr << "CXXLAPACK: WV_generic = " << WV_generic << std::endl;
//        std::cerr << "F77LAPACK: WV = " << WV << std::endl;
        failed = true;
    }

    if (! isIdentical(WH_generic, WH, "WH_generic", "WH")) {
//        std::cerr << "CXXLAPACK: WH_generic = " << WH_generic << std::endl;
//        std::cerr << "F77LAPACK: WH = " << WH << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr5.tcc" << std::endl;
        std::cerr << "N = H.numRows() =   " << H.numRows() << std::endl;
        std::cerr << "H.numCols() =       " << H.numCols() << std::endl;
        std::cerr << "NV = WV.numRows() = " << WV.numRows() << std::endl;
        std::cerr << "WV.numCols() =      " << WV.numCols() << std::endl;
        std::cerr << "NH = WH.numRows() = " << WH.numRows() << std::endl;
        std::cerr << "WH.numCols() =      " << WH.numCols() << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: laqr5.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename VSR, typename VSI, typename MH,
          typename MZ, typename MV, typename MU, typename MWV, typename MWH>
void
laqr5(bool                      wantT,
      bool                      wantZ,
      IndexType                 kacc22,
      IndexType                 kTop,
      IndexType                 kBot,
      IndexType                 nShifts,
      VSR                       &&sr,
      VSI                       &&si,
      MH                        &&H,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      MZ                        &&Z,
      MV                        &&V,
      MU                        &&U,
      MWV                       &&WV,
      MWH                       &&WH)
{
    CHECKPOINT_ENTER;
    laqr5(wantT, wantZ, kacc22, kTop, kBot, nShifts, sr, si, H, iLoZ, iHiZ, Z,
          V, U, WV, WH);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LAQR5_TCC
