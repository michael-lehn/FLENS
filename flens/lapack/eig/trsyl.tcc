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
       SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
      $                   LDC, SCALE, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_EIG_TRSYL_TCC
#define FLENS_LAPACK_EIG_TRSYL_TCC 1

#include <cmath>

namespace flens { namespace lapack {


//== generic lapack implementation =============================================
template <typename ISGN, typename MA, typename MB, typename MC, typename SCALE>
typename GeMatrix<MC>::IndexType
trsyl_generic(Transpose             transA,
              Transpose             transB,
              ISGN                  iSign,
              const GeMatrix<MA>    &A,
              const GeMatrix<MB>    &B,
              GeMatrix<MC>          &C,
              SCALE                 &scale)
{
    using std::abs;

    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    const T Zero(0), One(1);
//
//  .. Local Arrays ..
//
    T   _vecData[4], _xData[4];
    GeMatrixView<T>
        VEC    = typename GeMatrixView<T>::Engine(2, 2, _vecData, 2),
        X      = typename GeMatrixView<T>::Engine(2, 2, _xData, 2);

//
//  Decode and Test input parameters
//
    const bool noTransA = (transA==NoTrans);
    const bool noTransB = (transB==NoTrans);

    IndexType info = 0;
//
//  Quick return if possible
//
    scale = One;
    if (m==0 || n==0) {
        return info;
    }
//
//  Set constants to control overflow
//
    const T eps = lamch<T>(Precision);
    T smallNum = lamch<T>(SafeMin);
    T bigNum = One / smallNum;
    labad(smallNum, bigNum);
    smallNum = smallNum*T(m*n) / eps;
    bigNum = One / smallNum;

    const T sMin = max(smallNum,
                       eps*lan(MaximumNorm, A),
                       eps*lan(MaximumNorm, B));

    const T sign = iSign;

    IndexType   lNext, l1, l2, kNext, k1, k2;

    T           sumL, sumR, scaleC, A11, DA11, DB, xNorm;

    if (noTransA && noTransB) {
//
//      Solve    A*X + ISGN*X*B = scale*C.
//
//      The (K,L)th block of X is determined starting from
//      bottom-left corner column by column by
//
//       A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//
//      Where
//                 M                         L-1
//       R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
//               I=K+1                       J=1
//
//      Start column loop (index = L)
//      L1 (L2) : column index of the first (first) row of X(K,L).
//
        lNext = 1;
        for (IndexType l=1; l<=n; ++l) {
            if (l<lNext) {
                continue;
            }
            if (l==n) {
                l1 = l2 = l;
            } else {
                if (B(l+1,l)!=Zero) {
                    l1 = l;
                    l2 = l + 1;
                    lNext = l + 2;
                } else {
                    l1 = l2 = l;
                    lNext = l + 1;
                }
            }
//
//          Start row loop (index = K)
//          K1 (K2): row index of the first (last) row of X(K,L).
//
            kNext = m;
            for (IndexType k=m; k>=1; --k) {
                if (k>kNext) {
                    continue;
                }
                if (k==1) {
                    k1 = k2 = k;
                } else {
                    if (A(k,k-1)!=Zero) {
                        k1 = k - 1;
                        k2 = k;
                        kNext = k - 2;
                    } else {
                        k1 = k2 = k;
                        kNext = k - 1;
                    }
                }

                if (l1==l2 && k1==k2) {

                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);
                    scaleC = One;

                    A11 = A(k1,k1) + sign*B(l1,l1);
                    DA11 = abs(A11);
                    if (DA11<=sMin) {
                        A11 = sMin;
                        DA11 = sMin;
                        info = 1;
                    }
                    DB = abs(VEC(1,1));
                    if (DA11<One && DB>One) {
                        if (DB>bigNum*DA11 ) {
                            scaleC = One / DB;
                        }
                    }
                    X(1,1) = (VEC(1,1)*scaleC) / A11;

                    if (scaleC!=One) {
                        C     *=  scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);

                } else if (l1==l2 && k1!=k2) {

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    IndexType iErr = laln2(false, 1, sMin, One,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*B(l1,l1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *=  scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k2,l1) = X(2,1);

                } else if (l1!=l2 && k1==k2) {

                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = sign*(C(k1,l1)-(sumL+sign*sumR));

                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l2);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(2,1) = sign*(C(k1,l2)-(sumL+sign*sumR));

                    IndexType iErr = laln2(true, 1, sMin, One,
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*A(k1,k1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *=  scaleC;
                        scale *= scaleC;
                    }

                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(2,1);

                } else if (l1!=l2 && k1!=k2) {

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l2);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(1,2) = C(k1,l2) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l2);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(2,2) = C(k2,l2) - (sumL+sign*sumR);

                    IndexType iErr = lasy2(false, false, iSign,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           VEC, scaleC, X, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *=  scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(1,2);
                    C(k2,l1) = X(2,1);
                    C(k2,l2) = X(2,2);
                }

            }

        }

    } else if(!noTransA && noTransB) {
//
//      Solve    A**T *X + ISGN*X*B = scale*C.
//
//      The (K,L)th block of X is determined starting from
//      upper-left corner column by column by
//
//        A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//
//      Where
//                  K-1                          L-1
//         R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
//                  I=1                          J=1
//
//      Start column loop (index = L)
//      L1 (L2): column index of the first (last) row of X(K,L)
//
        lNext = 1;
        for (IndexType l=1; l<=n; ++l) {
            if (l<lNext) {
                continue;
            }
            if (l==n){
                l1 = l2 = l;
            } else {
                if (B(l+1,l)!=Zero) {
                    l1 = l;
                    l2 = l + 1;
                    lNext = l + 2;
                } else {
                    l1 = l2 = l;
                    lNext = l + 1;
                }
            }
//
//          Start row loop (index = K)
//          K1 (K2): row index of the first (last) row of X(K,L)
//
            kNext = 1;
            for (IndexType k=1; k<=m; ++k) {
                if (k<kNext) {
                    continue;
                }
                if (k==m) {
                    k1 = k2 = k;
                } else {
                    if (A(k+1,k)!=Zero) {
                        k1 = k;
                        k2 = k + 1;
                        kNext = k + 2;
                    } else {
                        k1 = k2 = k;
                        kNext = k + 1;
                    }
                }

                if (l1==l2 && k1==k2) {
                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);
                    scaleC = One;

                    A11 = A(k1,k1) + sign*B(l1,l1);
                    DA11 = abs(A11);
                    if (DA11<=sMin) {
                        A11 = sMin;
                        DA11 = sMin;
                        info = 1;
                    }
                    DB = abs(VEC(1,1));
                    if (DA11<One && DB>One) {
                        if (DB>bigNum*DA11) {
                            scaleC = One/DB;
                        }
                    }
                    X(1,1) = (VEC(1,1)*scaleC) / A11;

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);

                } else if (l1==l2 && k1!=k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l1);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    IndexType iErr = laln2(true, 1, sMin, One,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*B(l1,l1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k2,l1) = X(2,1);

                } else if (l1!=l2 && k1==k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = sign*(C(k1,l1)-(sumL+sign*sumR));

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l2);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(2,1) = sign*(C(k1,l2)-(sumL+sign*sumR));

                    IndexType iErr = laln2(true, 1, sMin, One,
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*A(k1,k1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(2,1);

                } else if (l1!=l2 && k1!=k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l2);
                    sumR = C(k1,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(1,2) = C(k1,l2) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l1);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l1);
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l2);
                    sumR = C(k2,_(1,l1-1)) * B(_(1,l1-1),l2);
                    VEC(2,2) = C(k2,l2) - (sumL+sign*sumR);

                    IndexType iErr = lasy2(true, false, iSign,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           VEC, scaleC, X, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(1,2);
                    C(k2,l1) = X(2,1);
                    C(k2,l2) = X(2,2);
                }

            }
        }

    } else if (!noTransA && !noTransB) {
//
//      Solve    A**T*X + ISGN*X*B**T = scale*C.
//
//      The (K,L)th block of X is determined starting from
//      top-right corner column by column by
//
//         A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
//
//      Where
//                   K-1                            N
//          R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
//                   I=1                          J=L+1
//
//      Start column loop (index = L)
//      L1 (L2): column index of the first (last) row of X(K,L)
//
        lNext = n;
        for (IndexType l=n; l>=1; --l) {
            if (l>lNext) {
                continue;
            }
            if (l==1) {
                l1 = l2 = l;
            } else {
                if (B(l,l-1)!=Zero) {
                    l1 = l-1;
                    l2 = l;
                    lNext = l - 2;
                } else {
                    l1 = l2 = l;
                    lNext = l - 1;
                }
            }
//
//          Start row loop (index = K)
//          K1 (K2): row index of the first (last) row of X(K,L)
//
            kNext = 1;
            for (IndexType k=1; k<=m; ++k) {
                if (k<kNext) {
                    continue;
                }
                if (k==m) {
                    k1 = k2 = k;
                } else {
                    if (A(k+1,k)!=Zero) {
                        k1 = k;
                        k2 = k + 1;
                        kNext = k + 2;
                    } else {
                        k1 = k2 = k;
                        kNext = k + 1;
                    }
                }

                if (l1==l2 && k1==k2) {
                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(l1+1,n)) * B(l1,_(l1+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);
                    scaleC = One;

                    A11 = A(k1,k1) + sign*B(l1,l1);
                    DA11 = abs(A11);
                    if (DA11<=sMin) {
                        A11 = sMin;
                        DA11 = sMin;
                        info = 1;
                    }
                    DB = abs(VEC(1,1));
                    if (DA11<One && DB>One) {
                        if (DB>bigNum*DA11) {
                            scaleC = One / DB;
                        }
                    }
                    X(1,1) = (VEC(1,1)*scaleC) / A11;

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);

                } else if (l1==l2 && k1!=k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l1);
                    sumR = C(k2,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    IndexType iErr = laln2(true, 1, sMin, One,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*B(l1,l1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k2,l1) = X(2,1);

                } else if (l1!=l2 && k1==k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = sign*(C(k1,l1)-(sumL+sign*sumR));

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l2);
                    sumR = C(k1,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(2,1) = sign*(C(k1,l2)-(sumL+sign*sumR));

                    IndexType iErr = laln2(false, 1, sMin, One,
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*A(k1,k1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(2,1);

                } else if (l1!=l2 && k1!=k2) {

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k1) * C(_(1,k1-1),l2);
                    sumR = C(k1,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(1,2) = C(k1,l2) - (sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l1);
                    sumR = C(k2,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(2,1) = C(k2,l1) - ( sumL+sign*sumR);

                    sumL = A(_(1,k1-1),k2) * C(_(1,k1-1),l2);
                    sumR = C(k2,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(2,2) = C(k2,l2) - (sumL+sign*sumR);

                    IndexType iErr = lasy2(true, true, iSign,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           VEC, scaleC, X, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(1,2);
                    C(k2,l1) = X(2,1);
                    C(k2,l2) = X(2,2);
                }

            }
        }

    } else if (noTransA && !noTransB) {
//
//      Solve    A*X + ISGN*X*B**T = scale*C.
//
//      The (K,L)th block of X is determined starting from
//      bottom-right corner column by column by
//
//          A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
//
//      Where
//                    M                          N
//          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
//                  I=K+1                      J=L+1
//
//      Start column loop (index = L)
//      L1 (L2): column index of the first (last) row of X(K,L)
//
        lNext = n;
        for (IndexType l=n; l>=1; --l) {
            if (l>lNext) {
                continue;
            }
            if (l==1) {
                l1 = l2 = l;
            } else {
                if (B(l,l-1)!=Zero) {
                    l1 = l - 1;
                    l2 = l;
                    lNext = l - 2;
                } else {
                    l1 = l2 = l;
                    lNext = l - 1;
                }
            }
//
//          Start row loop (index = K)
//          K1 (K2): row index of the first (last) row of X(K,L)
//
            kNext = m;
            for (IndexType k=m; k>=1; --k) {
                if (k>kNext) {
                    continue;
                }
                if (k==1) {
                    k1 = k2 = k;
                } else {
                    if (A(k,k-1)!=Zero) {
                        k1 = k - 1;
                        k2 = k;
                        kNext = k - 2;
                    } else {
                        k1 = k2 = k;
                        kNext = k - 1;
                    }
                }

                if (l1==l2 && k1==k2) {
                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l1);
                    sumR = C(k1,_(l1+1,n)) * B(l1,_(l1+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);
                    scaleC = One;

                    A11 = A(k1,k1) + sign*B(l1,l1);
                    DA11 = abs(A11);
                    if (DA11<=sMin) {
                        A11 = sMin;
                        DA11 = sMin;
                        info = 1;
                    }
                    DB = abs(VEC(1,1));
                    if (DA11<One && DB>One) {
                        if (DB>bigNum*DA11) {
                            scaleC = One / DB;
                        }
                    }
                    X(1,1) = (VEC(1,1)*scaleC) / A11;

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);

                } else if (l1==l2 && k1!=k2) {

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k2,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    IndexType iErr = laln2(false, 1, sMin, One,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*B(l1,l1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k2,l1) = X(2,1);

                } else if (l1!=l2 && k1==k2) {

                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = sign*(C(k1,l1)-(sumL+sign*sumR));

                    sumL = A(k1,_(k1+1,m)) * C(_(k1+1,m),l2);
                    sumR = C(k1,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(2,1) = sign*(C(k1,l2)-(sumL+sign*sumR));

                    IndexType iErr = laln2(false, 1, sMin, One,
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           One, One, VEC(_,_(1,1)),
                                           -sign*A(k1,k1), Zero,
                                           X(_,_(1,1)), scaleC, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(2,1);

                } else if (l1!=l2 && k1!=k2) {

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k1,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(1,1) = C(k1,l1) - (sumL+sign*sumR);

                    sumL = A(k1,_(k2+1,m)) * C(_(k2+1,m),l2);
                    sumR = C(k1,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(1,2) = C(k1,l2) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l1);
                    sumR = C(k2,_(l2+1,n)) * B(l1,_(l2+1,n));
                    VEC(2,1) = C(k2,l1) - (sumL+sign*sumR);

                    sumL = A(k2,_(k2+1,m)) * C(_(k2+1,m),l2);
                    sumR = C(k2,_(l2+1,n)) * B(l2,_(l2+1,n));
                    VEC(2,2) = C(k2,l2) - (sumL+sign*sumR);

                    IndexType iErr = lasy2(false, true, iSign,
                                           A(_(k1,k1+1),_(k1,k1+1)),
                                           B(_(l1,l1+1),_(l1,l1+1)),
                                           VEC, scaleC, X, xNorm);

                    if (iErr!=0) {
                        info = 1;
                    }

                    if (scaleC!=One) {
                        C     *= scaleC;
                        scale *= scaleC;
                    }
                    C(k1,l1) = X(1,1);
                    C(k1,l2) = X(1,2);
                    C(k2,l1) = X(2,1);
                    C(k2,l2) = X(2,2);
                }

            }
        }

    }
    return info;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename ISGN, typename MA, typename MB, typename MC, typename SCALE>
typename GeMatrix<MC>::IndexType
trsyl(Transpose             transA,
      Transpose             transB,
      ISGN                  sign,
      const GeMatrix<MA>    &A,
      const GeMatrix<MB>    &B,
      GeMatrix<MC>          &C,
      SCALE                 &scale)
{
    typedef typename GeMatrix<MC>::IndexType  IndexType;

    IndexType  info;
    info = cxxlapack::trsyl<IndexType>(getF77Char(transA),
                                       getF77Char(transB),
                                       sign,
                                       C.numRows(),
                                       C.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       B.data(),
                                       B.leadingDimension(),
                                       C.data(),
                                       C.leadingDimension(),
                                       scale);

    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename ISGN, typename MA, typename MB, typename MC, typename SCALE>
typename GeMatrix<MC>::IndexType
trsyl(Transpose             transA,
      Transpose             transB,
      ISGN                  iSign,
      const GeMatrix<MA>    &A,
      const GeMatrix<MB>    &B,
      GeMatrix<MC>          &C,
      SCALE                 &scale)
{
    LAPACK_DEBUG_OUT("trsyl");

    typedef typename GeMatrix<MC>::IndexType    IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(iSign==1 || iSign==-1);
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==B.numCols());

    const IndexType m = A.numRows();
    const IndexType n = B.numRows();

    ASSERT(C.firstRow()==1);
    ASSERT(C.firstCol()==1);
    ASSERT(C.numRows()==m);
    ASSERT(C.numCols()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    const typename GeMatrix<MC>::NoView   C_org       = C;
    SCALE                                 scale_org   = scale;
#   endif

//
//  Call implementation
//
    IndexType info = trsyl_generic(transA, transB, iSign, A, B, C, scale);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by generic implementation
//
    const typename GeMatrix<MC>::NoView   C_generic       = C;
    SCALE                                 scale_generic   = scale;
//
//  restore output arguments
//
    C       = C_org;
    scale   = scale_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::trsyl(transA, transB, iSign, A, B, C, scale);

    bool failed = false;
    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }
    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = "
                  << scale << std::endl;
        failed = true;
    }
    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }
    if (failed) {
        std::cerr << "transA = " << transA << std::endl;
        std::cerr << "transB = " << transB << std::endl;
        std::cerr << "iSign =  " << iSign << std::endl;
        std::cerr << "A =  " << A << std::endl;
        std::cerr << "B =  " << B << std::endl;

        std::cerr << "error in: trsyl.tcc" << std::endl;
        ASSERT(0);
    } else {
//      std::cerr << "passed: trsyl.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename ISGN, typename MA, typename MB, typename MC, typename SCALE>
typename MC::IndexType
trsyl(Transpose             transA,
      Transpose             transB,
      ISGN                  iSign,
      const GeMatrix<MA>    &A,
      const GeMatrix<MB>    &B,
      MC                    &&C,
      SCALE                 &&scale)
{
    typedef  typename GeMatrix<MC>::IndexType   IndexType;

    CHECKPOINT_ENTER;
    IndexType info = trsyl(transA, transB, iSign, A, B, C, scale);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_TRSYL_TCC
