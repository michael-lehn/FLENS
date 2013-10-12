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
      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, MM, M, WORK, INFO )
      SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, MM, M, WORK, RWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_TREVC_TCC
#define FLENS_LAPACK_IMPL_TREVC_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename IndexType, typename VWORK>
void
trevc_impl(bool                          computeVL,
           bool                          computeVR,
           TREVC::Job                    howMany,
           DenseVector<VSELECT>          &select,
           const GeMatrix<MT>            &T,
           GeMatrix<MVL>                 &VL,
           GeMatrix<MVR>                 &VR,
           IndexType                     mm,
           IndexType                     &m,
           DenseVector<VWORK>            &work)
{
    using std::abs;
    using flens::max;

    typedef typename GeMatrix<MT>::ElementType  ElementType;

    // TODO: define a colmajor GeMatrixView
    typedef typename GeMatrix<MT>::View         GeMatrixView;

    const ElementType   Zero(0), One(1);

    const Underscore<IndexType> _;
    const IndexType n = T.numCols();

//
//    .. Local Arrays ..
//
    ElementType  _xData[4];
    GeMatrixView X = typename GeMatrixView::Engine(2, 2, _xData, 2);

    GeMatrixView Work(n, 3, work, n);

//
//  Decode and test the input parameters
//
    const bool over     = howMany==TREVC::Backtransform;
    const bool someV    = howMany==TREVC::Selected;

//
//  Set M to the number of columns required to store the selected
//  eigenvectors, standardize the array SELECT if necessary, and
//  test MM.
//
    if (someV) {
        m = 0;
        bool pair = false;
        for (IndexType j=1; j<=n; ++j) {
            if (pair) {
                pair        = false;
                select(j)   = false;
            } else {
                if (j<n) {
                    if (T(j+1,j)==Zero) {
                        if (select(j)) {
                            ++m;
                        }
                    } else {
                        pair = true;
                        if (select(j) || select(j+1)) {
                            select(j) = true;
                            m += 2;
                        }
                    }
                } else {
                    if (select(n)) {
                        ++m;
                    }
                }
            }
        }
    } else {
        m = n;
    }

    if (mm<m) {
        ASSERT(0);
    }
//
//   Quick return if possible.
//
    if (n==0) {
        return;
    }
//
//  Set the constants to control overflow.
//
    ElementType underflow = lamch<ElementType>(SafeMin);
    ElementType overflow  = One / underflow;
    labad(underflow, overflow);

    const ElementType ulp         = lamch<ElementType>(Precision);
    const ElementType smallNum    = underflow*(n/ulp);
    const ElementType bigNum      = (One-ulp)/smallNum;
//
//  Compute 1-norm of each column of strictly upper triangular
//  part of T to control overflow in triangular solver.
//
    work(1) = Zero;
    for (IndexType j=2; j<=n; ++j) {
        work(j) = Zero;
        for (IndexType i=1; i<=j-1; ++i) {
            work(j) += abs(T(i,j));
        }
    }
//
//  Index IP is used to specify the real or complex eigenvalue:
//    IP = 0, real eigenvalue,
//         1, first of conjugate complex pair: (wr,wi)
//        -1, second of conjugate complex pair: (wr,wi)
//
    const IndexType n2 = 2*n;

    ElementType scale, xNorm, wr, wi;

    if (computeVR) {
//
//      Compute right eigenvectors.
//

        IndexType ip = 0;
        IndexType is = m;
        for (IndexType ki=n; ki>=1; --ki) {

            if (ip==1) {
                ip = 0;
                continue;
            }
            if (ki!=1 && T(ki,ki-1)!=Zero) {
                ip = -1;
            }

            if (someV) {
                if (ip==0) {
                    if (!select(ki)) {
                        if (ip==1) {
                            ip = 0;
                        }
                        if (ip==-1) {
                            ip = 1;
                        }
                        continue;
                    }
                } else {
                    if (!select(ki-1)) {
                        if (ip==1) {
                            ip = 0;
                        }
                        if (ip==-1) {
                            ip = 1;
                        }
                        continue;
                    }
                }
            }
//
//          Compute the KI-th eigenvalue (WR,WI).
//
            wr = T(ki,ki);
            wi = Zero;
            if (ip!=0) {
                wi = sqrt(abs(T(ki,ki-1)))*sqrt(abs(T(ki-1,ki)));
            }
            const ElementType safeMin = max(ulp*(abs(wr)+abs(wi)), smallNum);

            if (ip==0) {
//
//              Real right eigenvector
//
                work(ki+n) = One;
//
//              Form right-hand side
//              TODO: do this by work(_(...)) = -T(_(...),ki)
//
                for (IndexType k=1; k<=ki-1; ++k) {
                    work(k+n) = -T(k,ki);
                }
//
//              Solve the upper quasi-triangular system:
//                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
//
                IndexType jNext = ki - 1;
                for (IndexType j=ki-1; j>=1; --j) {
                    if (j>jNext) {
                        continue;
                    }
                    IndexType j1 = j;
                    IndexType j2 = j;
                    jNext = j - 1;
                    if (j>1) {
                        if (T(j,j-1)!=Zero) {
                            j1 = j - 1;
                            jNext = j - 2;
                        }
                    }

                    if (j1==j2) {
//
//                      1-by-1 diagonal block
//
                        laln2(false,
                              IndexType(1),
                              safeMin,
                              One,
                              T(_(j,j),_(j,j)),
                              One,
                              One,
                              Work(_(j,j), _(2,2)),
                              wr,
                              Zero,
                              X(_(1,1),_(1,1)),
                              scale,
                              xNorm);
//
//                      Scale X(1,1) to avoid overflow when updating
//                      the right-hand side.
//
                        if (xNorm>One) {
                            if (work(j)>bigNum/xNorm) {
                                X(1,1)  /= xNorm;
                                scale   /= xNorm;
                            }
                        }
//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(1+n,ki+n)) *= scale;
                        }
                        work(j+n) = X(1,1);
//
//                      Update right-hand side
//
                        work(_(1+n,j-1+n)) -= X(1,1)*T(_(1,j-1),j);

                    } else {
//
//                      2-by-2 diagonal block
//
                        laln2(false,
                              IndexType(1),
                              safeMin,
                              One,
                              T(_(j-1,j),_(j-1,j)),
                              One,
                              One,
                              Work(_(j-1,j), _(2,2)),
                              wr,
                              Zero,
                              X(_(1,2),_(1,1)),
                              scale,
                              xNorm);

//
//                      Scale X(1,1) and X(2,1) to avoid overflow when
//                      updating the right-hand side.
//
                        if (xNorm>One) {
                            const ElementType beta = max(work(j-1), work(j));
                            if (beta>bigNum/xNorm) {
                                X(1,1) /= xNorm;
                                X(2,1) /= xNorm;
                                scale  /= xNorm;
                            }
                        }
//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(1+n,ki+n)) *= scale;
                        }
                        work(j-1+n) = X(1,1);
                        work(j+n)   = X(2,1);
//
//                      Update right-hand side
//
                        work(_(1+n,j-2+n)) -= X(1,1)*T(_(1,j-2),j-1);
                        work(_(1+n,j-2+n)) -= X(2,1)*T(_(1,j-2),j);
                    }
                }
//
//              Copy the vector x or Q*x to VR and normalize.
//
                if (!over) {
                    VR(_(1,ki),is) = work(_(1+n,ki+n));

                    const IndexType ii = blas::iamax(VR(_(1,ki),is));
                    const ElementType reMax = One / abs(VR(ii,is));
                    VR(_(1,ki),is) *= reMax;

                    VR(_(ki+1,n),is) = Zero;
                } else {
                    if (ki>1) {
                        blas::mv(NoTrans, One,
                                 VR(_,_(1,ki-1)),
                                 work(_(1+n,n+ki-1)),
                                 work(ki+n),
                                 VR(_,ki));
                    }
                    const IndexType ii = blas::iamax(VR(_,ki));
                    const ElementType reMax = One / abs(VR(ii,ki));
                    VR(_,ki) *= reMax;
                }

            } else {
//
//              Complex right eigenvector.
//
//              Initial solve
//                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
//                [ (T(KI,KI-1)   T(KI,KI)   )               ]
//
                if (abs(T(ki-1,ki))>=abs(T(ki,ki-1))) {
                    work(ki-1+n)  = One;
                    work(ki+n2)   = wi / T(ki-1,ki);
                } else {
                    work(ki-1+n ) = -wi / T(ki,ki-1);
                    work(ki+n2 )  = One;
                }
                work(ki+n) = Zero;
                work(ki-1+n2) = Zero;
//
//              Form right-hand side
//
                for (IndexType k=1; k<=ki-2; ++k) {
                    work(k+n)   = -work(ki-1+n) * T(k,ki-1);
                    work(k+n2)  = -work(ki+n2)  * T(k,ki);
                }
//
//              Solve upper quasi-triangular system:
//              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
//
                IndexType jNext = ki - 2;
                for (IndexType j=ki-2; j>=1; --j) {
                    if (j>jNext) {
                        continue;
                    }
                    IndexType j1 = j;
                    IndexType j2 = j;
                    jNext = j - 1;
                    if (j>1) {
                        if (T(j,j-1)!=Zero) {
                            j1 = j - 1;
                            jNext = j - 2;
                        }
                    }

                    if (j1==j2) {
//
//                      1-by-1 diagonal block
//
                        laln2(false,
                              IndexType(2),
                              safeMin,
                              One,
                              T(_(j,j),_(j,j)),
                              One,
                              One,
                              Work(_(j,j), _(2,3)),
                              wr,
                              wi,
                              X(_(1,1),_(1,2)),
                              scale,
                              xNorm);
//
//                      Scale X(1,1) and X(1,2) to avoid overflow when
//                      updating the right-hand side.
//
                        if (xNorm>One) {
                            if (work(j)>bigNum/xNorm) {
                                X(1,1) /= xNorm;
                                X(1,2) /= xNorm;
                                scale /= xNorm;
                            }
                        }
//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(1+n,ki+n))   *= scale;
                            work(_(1+n2,ki+n2)) *= scale;
                        }
                        work(j+n)   = X(1,1);
                        work(j+n2)  = X(1,2);
//
//                      Update the right-hand side
//
                        work(_(1+n,j-1+n))   -= X(1,1)*T(_(1,j-1),j);
                        work(_(1+n2,j-1+n2)) -= X(1,2)*T(_(1,j-1),j);

                    } else {
//
//                      2-by-2 diagonal block
//
                        laln2(false,
                              IndexType(2),
                              safeMin,
                              One,
                              T(_(j-1,j),_(j-1,j)),
                              One,
                              One,
                              Work(_(j-1,j), _(2,3)),
                              wr,
                              wi,
                              X(_(1,2),_(1,2)),
                              scale,
                              xNorm);

//
//                      Scale X to avoid overflow when updating
//                      the right-hand side.
//
                        if (xNorm>One) {
                            const ElementType beta = max(work(j-1), work(j));
                            if (beta>bigNum/xNorm) {
                                const ElementType rec = One / xNorm;
                                X(1,1) *= rec;
                                X(1,2) *= rec;
                                X(2,1) *= rec;
                                X(2,2) *= rec;
                                scale *= rec;
                            }
                        }
//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(1+n,ki+n))   *= scale;
                            work(_(1+n2,ki+n2)) *= scale;
                        }
                        work(j-1+n)  = X(1,1);
                        work(j+n)    = X(2,1);
                        work(j-1+n2) = X(1,2);
                        work(j+n2)   = X(2,2);
//
//                      Update the right-hand side
//
                        work(_(1+n,j-2+n)) -= X(1,1)*T(_(1,j-2),j-1);
                        work(_(1+n,j-2+n)) -= X(2,1)*T(_(1,j-2),j);

                        work(_(1+n2,j-2+n2)) -= X(1,2)*T(_(1,j-2),j-1);
                        work(_(1+n2,j-2+n2)) -= X(2,2)*T(_(1,j-2),j);
                    }
                }
//
//              Copy the vector x or Q*x to VR and normalize.
//
                if (!over) {
                    VR(_(1,ki),is-1) = work(_(1+n,ki+n));
                    VR(_(1,ki),is)   = work(_(1+n2,ki+n2));

                    ElementType eMax = Zero;
                    for (IndexType k=1; k<=ki; ++k) {
                        eMax = max(eMax, abs(VR(k,is-1))+abs(VR(k,is)));
                    }

                    const ElementType reMax = One / eMax;
                    VR(_(1,ki),is-1)  *= reMax;
                    VR(_(1,ki),is)    *= reMax;

                    VR(_(ki+1,n),is-1)  = Zero;
                    VR(_(ki+1,n),is)    = Zero;

                } else {

                    if (ki>2) {
                        blas::mv(NoTrans, One,
                                 VR(_,_(1,ki-2)),
                                 work(_(1+n,ki-2+n)),
                                 work(ki-1+n),
                                 VR(_,ki-1));
                        blas::mv(NoTrans, One,
                                 VR(_,_(1,ki-2)),
                                 work(_(n2+1,n2+ki-2)),
                                 work(ki+n2),
                                 VR(_,ki));
                    } else {
                        VR(_,ki-1)  *= work(ki-1+n);
                        VR(_,ki)    *= work(ki+n2);
                    }

                    ElementType eMax = Zero;
                    for (IndexType k=1; k<=n; ++k) {
                        eMax = max(eMax, abs(VR(k,ki-1)) + abs(VR(k,ki)));
                    }
                    const ElementType reMax = One / eMax;
                    VR(_,ki-1)  *= reMax;
                    VR(_,ki)    *= reMax;
                }
            }

            --is;
            if (ip!=0) {
                --is;
            }
            if (ip==1) {
                ip = 0;
            }
            if (ip==-1) {
                ip = 1;
            }
        }
    }

    if (computeVL) {
//
//      Compute left eigenvectors.
//
        IndexType ip = 0;
        IndexType is = 1;
        for (IndexType ki=1; ki<=n; ++ki) {

            if (ip==-1) {
                ip = 0;
                continue;
            }
            if (ki!=n && T(ki+1,ki)!=Zero) {
                ip = 1;
            }

            if (someV) {
                if (!select(ki)) {
                    if (ip==-1) {
                        ip = 0;
                    }
                    if (ip==1) {
                        ip = -1;
                    }
                    continue;
                }
            }
//
//          Compute the KI-th eigenvalue (WR,WI).
//
            ElementType wr = T(ki,ki);
            ElementType wi = Zero;
            if (ip!=0) {
                wi = sqrt(abs(T(ki,ki+1))) * sqrt(abs(T(ki+1,ki)));
            }
            const ElementType safeMin = max(ulp*(abs(wr)+abs(wi)), smallNum);

            if (ip==0) {
//
//              Real left eigenvector.
//
                work(ki+n) = One;
//
//              Form right-hand side
//
                for (IndexType k=ki+1; k<=n; ++k) {
                    work(k+n) = -T(ki,k);
                }
//
//              Solve the quasi-triangular system:
//                 (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK
//
                ElementType vMax = One;
                ElementType vCrit = bigNum;

                IndexType jNext = ki + 1;
                for (IndexType j=ki+1; j<=n; ++j) {
                    if (j<jNext) {
                        continue;
                    }
                    IndexType j1 = j;
                    IndexType j2 = j;
                    jNext = j + 1;
                    if (j<n) {
                        if (T(j+1,j)!=Zero) {
                            j2 = j + 1;
                            jNext = j + 2;
                        }
                    }

                    if (j1==j2) {
//
//                      1-by-1 diagonal block
//
//                      Scale if necessary to avoid overflow when forming
//                      the right-hand side.
//
                        if (work(j)>vCrit) {
                            const ElementType rec = One / vMax;
                            work(_(ki+n,n2)) *= rec;
                            vMax = One;
                            vCrit = bigNum;
                        }

                        work(j+n) -= T(_(ki+1,j-1),j)*work(_(ki+1+n,n+j-1));
//
//                      Solve (T(J,J)-WR)**T*X = WORK
//
                        laln2(false,
                              IndexType(1),
                              safeMin,
                              One,
                              T(_(j,j),_(j,j)),
                              One,
                              One,
                              Work(_(j,j), _(2,2)),
                              wr,
                              Zero,
                              X(_(1,1),_(1,1)),
                              scale,
                              xNorm);

//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(ki+n,n2)) *= scale;
                        }
                        work(j+n) = X(1,1);
                        vMax = max(abs(work(j+n)), vMax);
                        vCrit = bigNum / vMax;

                    } else {
//
//                      2-by-2 diagonal block
//
//                      Scale if necessary to avoid overflow when forming
//                      the right-hand side.
//
                        const ElementType beta = max(work(j), work(j+1));
                        if (beta>vCrit) {
                            const ElementType rec = One / vMax;
                            work(_(ki+n,n2)) *= rec;
                            vMax = One;
                            vCrit = bigNum;
                        }

                        work(j+n)   -= T(_(ki+1,j-1),j)*work(_(ki+1+n,n+j-1));
                        work(j+1+n) -= T(_(ki+1,j-1),j+1)*work(_(ki+1+n,n+j-1));
//
//                      Solve
//                       [T(J,J)-WR   T(J,J+1)     ]**T * X = SCALE*( WORK1 )
//                       [T(J+1,J)    T(J+1,J+1)-WR]                ( WORK2 )
//
                        laln2(true,
                              IndexType(1),
                              safeMin,
                              One,
                              T(_(j,j+1),_(j,j+1)),
                              One,
                              One,
                              Work(_(j,j+1), _(2,2)),
                              wr,
                              Zero,
                              X(_(1,2),_(1,1)),
                              scale,
                              xNorm);

//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(ki+n,n2)) *= scale;
                        }
                        work(j+n)   = X(1,1);
                        work(j+1+n) = X(2,1);

                        vMax = max(abs(work(j+n)), abs(work(j+1+n)), vMax);
                        vCrit = bigNum / vMax;
                    }
                }
//
//              Copy the vector x or Q*x to VL and normalize.
//
                if (!over) {
                    VL(_(ki,n),is) = work(_(ki+n,n2));

                    const IndexType ii = blas::iamax(VL(_(ki,n),is)) + ki - 1;
                    const ElementType reMax = One / abs(VL(ii,is));
                    VL(_(ki,n),is) *= reMax;

                    VL(_(1,ki-1),is) = Zero;

                } else {

                    if (ki<n) {
                        blas::mv(NoTrans, One,
                                 VL(_,_(ki+1,n)),
                                 work(_(ki+1+n,n2)),
                                 work(ki+n),
                                 VL(_,ki));
                    }

                    const IndexType ii = blas::iamax(VL(_,ki));
                    const ElementType reMax = One / abs(VL(ii,ki));
                    VL(_,ki) *= reMax;

                }

            } else {
//
//              Complex left eigenvector.
//
//              Initial solve:
//                ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0.
//                ((T(KI+1,KI) T(KI+1,KI+1))                )
//
                if (abs(T(ki,ki+1))>=abs(T(ki+1,ki))) {
                    work(ki+n)    = wi / T(ki,ki+1);
                    work(ki+1+n2) = One;
                } else {
                    work(ki+n)    = One;
                    work(ki+1+n2) = -wi / T(ki+1,ki);
                }
                work(ki+1+n) = Zero;
                work(ki+n2)  = Zero;
//
//              Form right-hand side
//
                for (IndexType k=ki+2; k<=n; ++k) {
                    work(k+n)   = -work(ki+n)*T(ki,k);
                    work(k+n2)  = -work(ki+1+n2)*T(ki+1,k);
                }
//
//              Solve complex quasi-triangular system:
//              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
//
                ElementType vMax = One;
                ElementType vCrit = bigNum;

                IndexType jNext = ki + 2;
                for (IndexType j=ki+2; j<=n; ++j) {
                    if (j<jNext) {
                        continue;
                    }
                    IndexType j1 = j;
                    IndexType j2 = j;
                    jNext = j + 1;
                    if (j<n) {
                        if (T(j+1,j)!=Zero) {
                            j2 = j + 1;
                            jNext = j + 2;
                        }
                    }

                    if (j1==j2) {
//
//                      1-by-1 diagonal block
//
//                      Scale if necessary to avoid overflow when
//                      forming the right-hand side elements.
//
                        if (work(j)>vCrit) {
                            const ElementType rec = One / vMax;
                            work(_(ki+n,n2))    *= rec;
                            work(_(ki+n2,n2+n)) *= rec;
                            vMax = One;
                            vCrit = bigNum;
                        }

                        work(j+n)  -= T(_(ki+2,j-1),j)*work(_(ki+2+n,n+j-1));
                        work(j+n2) -= T(_(ki+2,j-1),j)*work(_(ki+2+n2,n2+j-1));
//
//                      Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
//
                        laln2(true,
                              IndexType(2),
                              safeMin,
                              One,
                              T(_(j,j),_(j,j)),
                              One,
                              One,
                              Work(_(j,j), _(2,3)),
                              wr,
                              -wi,
                              X(_(1,1),_(1,2)),
                              scale,
                              xNorm);

//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(ki+n,n2))    *= scale;
                            work(_(ki+n2,n2+n)) *= scale;
                        }
                        work(j+n) =  X(1,1);
                        work(j+n2) = X(1,2);
                        vMax = max(abs(work(j+n)), abs(work(j+n2)), vMax);
                        vCrit = bigNum / vMax;

                    } else {
//
//                      2-by-2 diagonal block
//
//                      Scale if necessary to avoid overflow when forming
//                      the right-hand side elements.
//
                        const ElementType beta = max(work(j), work(j+1));
                        if (beta>vCrit) {
                            const ElementType rec = One / vMax;
                            work(_(ki+n,n2))    *= rec;
                            work(_(ki+n2,n2+n)) *= rec;
                            vMax = One;
                            vCrit = bigNum;
                        }

                        const auto _work1 = work(_(ki+2+n,n+j-1));
                        const auto _work2 = work(_(ki+2+n2,n2+j-1));

                        work(j+n)   -= T(_(ki+2,j-1),j)   * _work1;
                        work(j+n2)  -= T(_(ki+2,j-1),j)   * _work2;
                        work(j+1+n) -= T(_(ki+2,j-1),j+1) * _work1;
                        work(j+1+n2)-= T(_(ki+2,j-1),j+1) * _work2;
//
//                      Solve 2-by-2 complex linear equation
//                      ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B
//                      ([T(j+1,j) T(j+1,j+1)]               )
//
                        laln2(true,
                              IndexType(2),
                              safeMin,
                              One,
                              T(_(j,j+1),_(j,j+1)),
                              One,
                              One,
                              Work(_(j,j+1), _(2,3)),
                              wr,
                              -wi,
                              X(_(1,2),_(1,2)),
                              scale,
                              xNorm);
//
//                      Scale if necessary
//
                        if (scale!=One) {
                            work(_(ki+n,n2))    *= scale;
                            work(_(ki+n2,n2+n)) *= scale;
                        }
                        work(j+n)    = X(1,1);
                        work(j+n2)   = X(1,2);
                        work(j+1+n)  = X(2,1);
                        work(j+1+n2) = X(2,2);
                        vMax = max(abs(X(1,1)), abs(X(1,2)),
                                   abs(X(2,1)), abs(X(2,2)),
                                   vMax);
                        vCrit = bigNum / vMax;

                    }
                }
//
//              Copy the vector x or Q*x to VL and normalize.
//
                if (!over) {
                    VL(_(ki,n),is)   = work(_(ki+n,n2));
                    VL(_(ki,n),is+1) = work(_(ki+n2,n2+n));

                    ElementType eMax = Zero;
                    for (IndexType k=ki; k<=n; ++k) {
                        eMax = max(eMax, abs(VL(k,is)) + abs(VL(k,is+1)));
                    }
                    const ElementType reMax = One / eMax;
                    VL(_(ki,n),is)   *= reMax;
                    VL(_(ki,n),is+1) *= reMax;

                    VL(_(1,ki-1),is)   = Zero;
                    VL(_(1,ki-1),is+1) = Zero;
                } else {
                    if (ki<n-1) {
                        blas::mv(NoTrans, One,
                                 VL(_,_(ki+2,n)),
                                 work(_(ki+2+n,n2)),
                                 work(ki+n),
                                 VL(_,ki));
                        blas::mv(NoTrans, One,
                                 VL(_,_(ki+2,n)),
                                 work(_(ki+2+n2,n2+n)),
                                 work(ki+1+n2),
                                 VL(_,ki+1));
                    } else {
                        VL(_,ki)    *= work(ki+n);
                        VL(_,ki+1)  *= work(ki+1+n2);
                    }

                    ElementType eMax = Zero;
                    for (IndexType k=1; k<=n; ++k) {
                        eMax = max(eMax, abs(VL(k,ki)) + abs(VL(k,ki+1)));
                    }
                    const ElementType reMax = One / eMax;
                    VL(_,ki)   *= reMax;
                    VL(_,ki+1) *= reMax;
                }
            }

            ++is;
            if (ip!=0) {
                ++is;
            }
            if (ip==-1) {
                ip = 0;
            }
            if (ip==1) {
                ip = -1;
            }
        }
    }
}

//
//  Complex variant
//
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename IndexType, typename VWORK, typename VRWORK>
void
trevc_impl(bool                          computeVL,
           bool                          computeVR,
           TREVC::Job                    howMany,
           DenseVector<VSELECT>          &select,
           GeMatrix<MT>                  &T,
           GeMatrix<MVL>                 &VL,
           GeMatrix<MVR>                 &VR,
           IndexType                     mm,
           IndexType                     &m,
           DenseVector<VWORK>            &work,
           DenseVector<VRWORK>           &rWork)
{
    using std::abs;
    using cxxblas::abs1;
    using flens::max;

    typedef typename GeMatrix<MT>::ElementType                  ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType   PrimitiveType;

    const PrimitiveType  Zero(0), One(1);
    const ElementType    CZero(0), COne(1);

    const Underscore<IndexType> _;
    const IndexType n = T.numCols();

//
//  Decode and test the input parameters
//
    const bool over     = howMany==TREVC::Backtransform;
    const bool someV    = howMany==TREVC::Selected;

//
//  Set M to the number of columns required to store the selected
//  eigenvectors.
//
    if (someV) {
        m = 0;
        for (IndexType j=1; j<=n; ++j) {
            if (select(j)) {
                ++m;
            }
        }
    } else {
        m = n;
    }

    if (mm<m) {
        ASSERT(0);
    }
//
//  Quick return if possible.
//
    if (n==0) {
        return;
    }
//
//  Set the constants to control overflow.
//
    PrimitiveType underflow = lamch<PrimitiveType>(SafeMin);
    PrimitiveType overflow  = One / underflow;
    labad(underflow, overflow);

    const PrimitiveType ulp         = lamch<PrimitiveType>(Precision);
    const PrimitiveType smallNum    = underflow*(n/ulp);
//
//  Store the diagonal elements of T in working array WORK.
//
    for (IndexType i=1; i<=n; ++i) {
        work(i+n) = T(i,i);
    }
//
//  Compute 1-norm of each column of strictly upper triangular
//  part of T to control overflow in triangular solver.
//
    rWork(1) = Zero;
    for (IndexType j=2; j<=n; ++j) {
        rWork(j) = blas::asum(T(_(1,j-1),j));
    }

    if (computeVR) {
//
//      Compute right eigenvectors.
//
        IndexType       is = m;
        PrimitiveType   scale;

        for (IndexType ki=n; ki>=1; --ki) {

            if (someV) {
                if (! select(ki)) {
                    break;
                }
            }
            const PrimitiveType safeMin = max(ulp*abs1(T(ki,ki)), smallNum);

            work(1) = COne;
//
//          Form right-hand side.
//
            for (IndexType k=1; k<=ki-1; ++k) {
                work(k) = -T(k,ki);
            }
//
//          Solve the triangular system:
//             (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
//
            for (IndexType k=1; k<=ki-1; ++k) {
                T(k,k) -= T(ki,ki);
                if (abs1(T(k,k)) < safeMin) {
                    T(k,k) = safeMin;
                }
            }

            if (ki>1) {
                const auto _T    = T(_(1,ki-1),_(1,ki-1));
                auto       _x    = work(_(1,ki-1));
                auto       _norm = rWork(_(1,ki-1));

                latrs(NoTrans, true, _T.upper(), _x, scale, _norm);
                work(ki) = scale;
            }
//
//          Copy the vector x or Q*x to VR and normalize.
//
            if (! over) {
                VR(_(1,ki),is) = work(_(1,ki));

                IndexType     ii   = blas::iamax(VR(_(1,ki),is));
                PrimitiveType rMax = One / abs1(VR(ii,is));

                blas::scal(rMax, VR(_(1,ki),is));

                for (IndexType k=ki+1; k<=n; ++k) {
                    VR(k,is) = CZero;
                }
            } else {
                if (ki>1) {
                    const ElementType cScale = scale;
                    const auto        _VR    = VR(_,_(1,ki-1));
                    const auto        _x     = work(_(1,ki-1));
                    blas::mv(NoTrans, COne, _VR, _x, cScale, VR(_,ki));
                }
                IndexType     ii   = blas::iamax(VR(_,ki));
                PrimitiveType rMax = One / abs1(VR(ii,ki));
                blas::scal(rMax, VR(_,ki));
            }
//
//          Set back the original diagonal elements of T.
//
            for (IndexType k=1; k<=ki-1; ++k) {
                T(k,k) = work(k+n);
            }

            --is;
        }
    }

    if (computeVL) {
//
//      Compute left eigenvectors.
//
        IndexType       is = 1;
        PrimitiveType   scale;

        for (IndexType ki=1; ki<=n; ++ki) {

            if (someV) {
                if (! select(ki)) {
                    break;
                }
            }
            const PrimitiveType safeMin = max(ulp*abs1(T(ki,ki)), smallNum);

            work(n) = COne;
//
//          Form right-hand side.
//
            for (IndexType k=ki+1; k<=n; ++k) {
                work(k) = -conj(T(ki,k));
            }
//
//          Solve the triangular system:
//             (T(KI+1:N,KI+1:N) - T(KI,KI))**H * X = SCALE*WORK.
//
            for (IndexType k=ki+1; k<=n; ++k) {
                T(k,k) -= T(ki,ki);
                if (abs1(T(k,k)) < safeMin) {
                    T(k,k) = safeMin;
                }
            }

            if (ki < n) {
                const auto _T    = T(_(ki+1,n),_(ki+1,n));
                auto       _x    = work(_(ki+1,n));
                auto       _norm = rWork(_(1,n-ki));

                latrs(ConjTrans, true, _T.upper(), _x, scale, _norm);
                work(ki) = scale;
            }
//
//          Copy the vector x or Q*x to VL and normalize.
//
            if (! over) {
                VL(_(ki,n),is) = work(_(ki,n));

                IndexType     ii   = blas::iamax(VL(_(ki,n),is)) + ki - 1;
                PrimitiveType rMax = One / abs1(VL(ii,is));

                blas::scal(rMax, VL(_(ki,n),is));

                for (IndexType k=1; k<=ki-1; ++k) {
                    VL(k,is) = CZero;
                }
            } else {
                if (ki<n) {
                    const ElementType cScale = scale;
                    const auto        _VL    = VL(_,_(ki+1,n));
                    const auto        _x     = work(_(ki+1,n));
                    blas::mv(NoTrans, COne, _VL, _x, cScale, VL(_,ki));
                }

                IndexType     ii = blas::iamax(VL(_,ki));
                PrimitiveType rMax = One / abs1(VL(ii,ki));
                blas::scal(rMax, VL(_,ki));
            }
//
//          Set back the original diagonal elements of T.
//
            for (IndexType k=ki+1; k<=n; ++k) {
                T(k,k) = work(k+n);
            }

            ++is;
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename IndexType, typename VWORK>
void
trevc_impl(bool                           computeVL,
           bool                           computeVR,
           TREVC::Job                     howMany,
           DenseVector<VSELECT>           &select,
           const GeMatrix<MT>             &T,
           GeMatrix<MVL>                  &VL,
           GeMatrix<MVR>                  &VR,
           IndexType                      mm,
           IndexType                      &m,
           DenseVector<VWORK>             &work)
{
    char side = 'N';
    if (computeVL && computeVR) {
        side = 'B';
    } else if (computeVL) {
        side = 'L';
    } else if (computeVR) {
        side = 'R';
    } else {
        ASSERT(0);
    }

    DenseVector<Array<IndexType> > _select = select;

    cxxlapack::trevc<IndexType>(side,
                                getF77Char(howMany),
                                _select.data(),
                                T.numRows(),
                                T.data(),
                                T.leadingDimension(),
                                VL.data(),
                                VL.leadingDimension(),
                                VR.data(),
                                VR.leadingDimension(),
                                mm,
                                m,
                                work.data());
}

//
//  Complex variant
//
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename IndexType, typename VWORK, typename VRWORK>
void
trevc_impl(bool                           computeVL,
           bool                           computeVR,
           TREVC::Job                     howMany,
           DenseVector<VSELECT>           &select,
           GeMatrix<MT>                   &T,
           GeMatrix<MVL>                  &VL,
           GeMatrix<MVR>                  &VR,
           IndexType                      mm,
           IndexType                      &m,
           DenseVector<VWORK>             &work,
           DenseVector<VRWORK>            &rWork)
{
    char side = 'N';
    if (computeVL && computeVR) {
        side = 'B';
    } else if (computeVL) {
        side = 'L';
    } else if (computeVR) {
        side = 'R';
    } else {
        ASSERT(0);
    }

    DenseVector<Array<IndexType> > _select = select;

    cxxlapack::trevc<IndexType>(side,
                                getF77Char(howMany),
                                _select.data(),
                                T.numRows(),
                                T.data(),
                                T.leadingDimension(),
                                VL.data(),
                                VL.leadingDimension(),
                                VR.data(),
                                VR.leadingDimension(),
                                mm,
                                m,
                                work.data(),
                                rWork.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//
//  Real variant
//
template <typename VSELECT, typename MTR, typename MVL, typename MVR,
          typename IndexType, typename VWORK>
typename RestrictTo<IsIntegerDenseVector<VSELECT>::value
                 && IsRealGeMatrix<MTR>::value
                 && IsRealGeMatrix<MVL>::value
                 && IsRealGeMatrix<MVR>::value
                 && IsInteger<IndexType>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
trevc(bool                          computeVL,
      bool                          computeVR,
      TREVC::Job                    howMany,
      VSELECT                       &select,
      const MTR                     &T,
      MVL                           &&VL,
      MVR                           &&VR,
      IndexType                     mm,
      IndexType                     &m,
      VWORK                         &&work)
{
    LAPACK_DEBUG_OUT("BEGIN: trevc");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<VSELECT>::Type   VectorSelect;
    typedef typename RemoveRef<MVL>::Type       MatrixVL;
    typedef typename RemoveRef<MVR>::Type       MatrixVR;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType n = T.numRows();

    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());
    ASSERT(VL.firstRow()==1);
    ASSERT(VL.firstCol()==1);
    ASSERT(VR.firstRow()==1);
    ASSERT(VR.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==3*n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename VectorSelect::NoView   select_org = select;
    typename MatrixVL::NoView       VL_org     = VL;
    typename MatrixVR::NoView       VR_org     = VR;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::trevc_impl(computeVL, computeVR, howMany, select,
                              T, VL, VR, mm, m, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename VectorSelect::NoView   select_generic = select;
    typename MatrixVL::NoView       VL_generic     = VL;
    typename MatrixVR::NoView       VR_generic     = VR;
    typename VectorWork::NoView     work_generic   = work;
//
//  restore output arguments
//
    select = select_org;
    VL     = VL_org;
    VR     = VR_org;
    work   = work_org;
//
//  Compare results
//
    external::trevc_impl(computeVL, computeVR, howMany, select,
                         T, VL, VR, mm, m, work);

    bool failed = false;
    if (! isIdentical(VL_generic, VL, "VL_generic", "VL")) {
        std::cerr << "CXXLAPACK: VL_generic = " << VL_generic << std::endl;
        std::cerr << "F77LAPACK: VL = " << VL << std::endl;
        failed = true;
    }

    if (! isIdentical(VR_generic, VR, "VR_generic", "_VR")) {
        std::cerr << "CXXLAPACK: VR_generic = " << VR_generic << std::endl;
        std::cerr << "F77LAPACK: VR = " << VR << std::endl;
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
    LAPACK_DEBUG_OUT("END: trevc");
}

//
//  Complex variant
//
template <typename VSELECT, typename MTR, typename MVL, typename MVR,
          typename IndexType, typename VWORK, typename VRWORK>
typename RestrictTo<IsIntegerDenseVector<VSELECT>::value
                 && IsComplexGeMatrix<MTR>::value
                 && IsComplexGeMatrix<MVL>::value
                 && IsComplexGeMatrix<MVR>::value
                 && IsInteger<IndexType>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         void>::Type
trevc(bool                          computeVL,
      bool                          computeVR,
      TREVC::Job                    howMany,
      VSELECT                       &select,
      MTR                           &&Tr,
      MVL                           &&VL,
      MVR                           &&VR,
      IndexType                     mm,
      IndexType                     &m,
      VWORK                         &&work,
      VRWORK                        &&rWork)
{
    LAPACK_DEBUG_OUT("BEGIN: trevc (complex)");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<VSELECT>::Type   VectorSelect;
    typedef typename RemoveRef<MTR>::Type       MatrixTr;
    typedef typename RemoveRef<MVL>::Type       MatrixVL;
    typedef typename RemoveRef<MVR>::Type       MatrixVR;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
    typedef typename RemoveRef<VRWORK>::Type    VectorRWork;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType n = Tr.numRows();

    ASSERT(Tr.firstRow()==1);
    ASSERT(Tr.firstCol()==1);
    ASSERT(Tr.numRows()==Tr.numCols());
    ASSERT(VL.firstRow()==1);
    ASSERT(VL.firstCol()==1);
    ASSERT(VR.firstRow()==1);
    ASSERT(VR.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==2*n);
    ASSERT(rWork.firstIndex()==1);
    ASSERT(rWork.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename VectorSelect::NoView   select_org = select;
    typename MatrixTr::NoView       Tr_org     = Tr;
    typename MatrixVL::NoView       VL_org     = VL;
    typename MatrixVR::NoView       VR_org     = VR;
    typename VectorWork::NoView     work_org   = work;
    typename VectorRWork::NoView    rWork_org  = rWork;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::trevc_impl(computeVL, computeVR, howMany, select,
                              Tr, VL, VR, mm, m, work, rWork);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename VectorSelect::NoView   select_generic = select;
    typename MatrixTr::NoView       Tr_generic     = Tr;
    typename MatrixVL::NoView       VL_generic     = VL;
    typename MatrixVR::NoView       VR_generic     = VR;
    typename VectorWork::NoView     work_generic   = work;
    typename VectorRWork::NoView    rWork_generic  = rWork;
//
//  restore output arguments
//
    select = select_org;
    Tr     = Tr_org;
    VL     = VL_org;
    VR     = VR_org;
    work   = work_org;
    rWork  = rWork_org;
//
//  Compare results
//
    external::trevc_impl(computeVL, computeVR, howMany, select,
                         Tr, VL, VR, mm, m, work, rWork);

    bool failed = false;
    if (! isIdentical(Tr_generic, Tr, "Tr_generic", "Tr")) {
        std::cerr << "CXXLAPACK: Tr_generic = " << Tr_generic << std::endl;
        std::cerr << "F77LAPACK: Tr = " << Tr << std::endl;
        failed = true;
    }

    if (! isIdentical(VL_generic, VL, "VL_generic", "_VL")) {
        std::cerr << "CXXLAPACK: VL_generic = " << VL_generic << std::endl;
        std::cerr << "F77LAPACK: VL = " << VL << std::endl;
        failed = true;
    }

    if (! isIdentical(VR_generic, VR, "VR_generic", "_VR")) {
        std::cerr << "CXXLAPACK: VR_generic = " << VR_generic << std::endl;
        std::cerr << "F77LAPACK: VR = " << VR << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(rWork_generic, rWork, "rWork_generic", "rWork")) {
        std::cerr << "CXXLAPACK: rWork_generic = " << rWork_generic
                  << std::endl;
        std::cerr << "F77LAPACK: rWork = " << rWork << std::endl;
        failed = true;
    }
    if (failed) {
        ASSERT(0);
    }
#   endif
    LAPACK_DEBUG_OUT("END: trevc (complex)");
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TREVC_TCC
