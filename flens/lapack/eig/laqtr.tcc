/*
 *  Copyright (c) 2011, Michael Lehn
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *  1) Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2) Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *  3) Neither the name of the FLENS development group nor the names of
 *     its contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Based on
      SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK,
     $                   INFO )
 *
 * -- LAPACK auxiliary routine (version 3.3.1) --
 * -- LAPACK is a software package provided by Univ. of Tennessee,    --
 * -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 * -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_EIG_LAQTR_TCC
#define FLENS_LAPACK_EIG_LAQTR_TCC 1

#include <flens/aux/aux.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MT, typename VB, typename W, typename SCALE, typename VX,
          typename VWORK>
typename DenseVector<VX>::IndexType
laqtr_generic(bool                  trans,
              bool                  real,
              const GeMatrix<MT>    &T,
              const DenseVector<VB> &b,
              const W               &w,
              SCALE                 &scale,
              DenseVector<VX>       &x,
              DenseVector<VWORK>    &work)
{
    using std::abs;
    using flens::max;

    typedef typename DenseVector<VX>::ElementType  ElementType;
    typedef typename DenseVector<VX>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const IndexType n = T.numRows();

    const ElementType Zero(0), One(1);
//
//  .. Local Arrays ..
//
    ElementType     _dataD[4], _dataV[4];
    GeMatrixView<ElementType>
        D = typename GeMatrixView<ElementType>::Engine(2, 2, _dataD, 2),
        V = typename GeMatrixView<ElementType>::Engine(2, 2, _dataV, 2);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Set constants to control overflow
//
    const ElementType eps = lamch<ElementType>(Precision);
    const ElementType smallNum = lamch<ElementType>(SafeMin)/eps;
    const ElementType bigNum = One/smallNum;

    ElementType xNorm = lan(MaximumNorm, T);
    if (!real) {
        typedef typename DenseVector<VB>::ElementType TB;
        const GeMatrixConstView<TB>  B(n, 1, b, n);

        xNorm = max(xNorm, abs(w), lan(MaximumNorm, B));
    }
    const ElementType sMin = max(smallNum, eps*xNorm);
//
//  Compute 1-norm of each column of strictly upper triangular
//  part of T to control overflow in triangular solver.
//
    work(1) = Zero;
    for (IndexType j=2; j<=n; ++j) {
        work(j) = blas::asum(T(_(1,j-1),j));
    }

    if (!real) {
        for (IndexType i=2; i<=n; ++i) {
            work(i) += abs(b(i));
        }
    }

    const IndexType n2 = 2*n;
    const IndexType n1 = (real) ? n : n2;

    const IndexType k = blas::iamax(x(_(1,n1)));
    ElementType xMax  = abs(x(k));
    scale = One;

    if (xMax>bigNum) {
        scale = bigNum / xMax;
        x(_(1,n1)) *= scale;
        xMax = bigNum;
    }

    if (real) {

        if (!trans) {
//
//          Solve T*p = scale*c
//
            IndexType jNext = n;
            for (IndexType j=n; j>=1; --j) {
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
//                  Meet 1 by 1 diagonal block
//
//                  Scale to avoid overflow when computing
//                      x(j) = b(j)/T(j,j)
//
                    ElementType xj  = abs(x(j1));
                    ElementType Tjj = abs(T(j1,j1));
                    ElementType tmp = T(j1,j1);
                    if (Tjj<sMin) {
                        tmp  = sMin;
                        Tjj  = sMin;
                        info = 1;
                    }

                    if (xj==Zero) {
                        continue;
                    }

                    if (Tjj<One) {
                        if (xj>bigNum*Tjj) {
                            const ElementType rec = One / xj;
                            x(_(1,n)) *= rec;
                            scale     *= rec;
                            xMax      *= rec;
                        }
                    }
                    x(j1) /= tmp;
                    xj = abs(x(j1));
//
//                  Scale x if necessary to avoid overflow when adding a
//                  multiple of column j1 of T.
//
                    if (xj>One) {
                        const ElementType rec = One / xj;
                        if (work(j1)>(bigNum-xMax)*rec) {
                            x(_(1,n)) *= rec;
                            scale     *= rec;
                        }
                    }
                    if (j1>1) {
                        auto _x = x(_(1,j1-1));

                        _x -= x(j1)*T(_(1,j1-1),j1);
                        const IndexType k = blas::iamax(_x);
                        xMax = abs(x(k));
                    }

                } else {
//
//                  Meet 2 by 2 diagonal block
//
//                  Call 2 by 2 linear system solve, to take
//                  care of possible overflow by scaling factor.
//
                    D(1,1) = x(j1);
                    D(2,1) = x(j2);

                    ElementType scaleC;
                    IndexType iErr = laln2(false, IndexType(1), sMin, One,
                                           T(_(j1,j1+1),_(j1,j1+1)), One, One,
                                           D(_(1,2),_(1,1)), Zero, Zero,
                                           V(_(1,2),_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 2;
                    }

                    if (scaleC!=One) {
                        x(_(1,n)) *= scaleC;
                        scale     *= scaleC;
                    }
                    x(j1) = V(1,1);
                    x(j2) = V(2,1);
//
//                  Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
//                  to avoid overflow in updating right-hand side.
//
                    const ElementType xj = max(abs(V(1,1)), abs(V(2,1)));
                    if (xj>One) {
                        const ElementType rec = One / xj;
                        if (max(work(j1),work(j2))>(bigNum-xMax)*rec) {
                            x(_(1,n)) *= rec;
                            scale     *= rec;
                        }
                    }
//
//                  Update right-hand side
//
                    if (j1>1) {
                        auto _x = x(_(1,j1-1));
                        _x -= x(j1) * T(_(1,j1-1),j1);
                        _x -= x(j2) * T(_(1,j1-1),j2);
                        const IndexType k = blas::iamax(_x);
                        xMax = abs(x(k));
                    }

                }

            }

        } else {
//
//          Solve T**T*p = scale*c
//
            IndexType jNext = 1;
            for (IndexType j=1; j<=n; ++j) {
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
//                  1 by 1 diagonal block
//
//                  Scale if necessary to avoid overflow in forming the
//                  right-hand side element by inner product.
//
                    ElementType xj = abs(x(j1));
                    if (xMax>One) {
                        const ElementType rec = One / xMax;
                        if (work(j1)>(bigNum-xj)*rec) {
                            x(_(1,n)) *= rec;
                            scale     *= rec;
                            xMax      *= rec;
                        }
                    }

                    x(j1) -= T(_(1,j1-1),j1) * x(_(1,j1-1));

                    xj = abs(x(j1));
                    ElementType Tjj = abs(T(j1,j1));
                    ElementType tmp = T(j1,j1);
                    if (Tjj<sMin) {
                        tmp = sMin;
                        Tjj = sMin;
                        info = 1;
                    }

                    if (Tjj<One) {
                        if (xj>bigNum*Tjj) {
                            const ElementType rec = One / xj;
                            x(_(1,n)) *= rec;
                            scale     *= rec;
                            xMax      *= rec;
                        }
                    }
                    x(j1) /= tmp;
                    xMax = max(xMax, abs(x(j1)));

                } else {
//
//                  2 by 2 diagonal block
//
//                  Scale if necessary to avoid overflow in forming the
//                  right-hand side elements by inner product.
//
                    const ElementType xj = max(abs(x(j1)), abs(x(j2)));
                    if (xMax>One) {
                        const ElementType rec = One / xMax;
                        if (max(work(j2),work(j1))>(bigNum-xj)*rec) {
                            x(_(1,n)) *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                    }

                    D(1,1) = x(j1) - T(_(1,j1-1),j1) * x(_(1,j1-1));
                    D(2,1) = x(j2) - T(_(1,j1-1),j2) * x(_(1,j1-1));

                    ElementType scaleC;
                    IndexType iErr = laln2(true, IndexType(1), sMin, One,
                                           T(_(j1,j1+1),_(j1,j1+1)), One, One,
                                           D(_(1,2),_(1,1)), Zero, Zero,
                                           V(_(1,2),_(1,1)), scaleC, xNorm);
                    if (iErr!=0) {
                        info = 2;
                    }

                    if (scaleC!=One) {
                        x(_(1,n)) *= scaleC;
                        scale     *= scaleC;
                    }
                    x(j1) = V(1,1);
                    x(j2) = V(2,1);
                    xMax = max(abs(x(j1)), abs(x(j2)), xMax);

                }
            }
        }

    } else {

        const ElementType sMinW = max(eps*abs(w), sMin);
        if (!trans) {
//
//          Solve (T + iB)*(p+iq) = c+id
//
            IndexType jNext = n;
            for (IndexType j=n; j>=1; --j) {
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
//                  1 by 1 diagonal block
//
//                  Scale if necessary to avoid overflow in division
//
                    ElementType z = w;
                    if (j1==1) {
                        z = b(1);
                    }
                    ElementType xj = abs(x(j1)) + abs(x(n+j1));
                    ElementType Tjj = abs(T(j1,j1)) + abs(z);
                    ElementType tmp = T(j1,j1);
                    if (Tjj<sMinW) {
                        tmp = sMinW;
                        Tjj = sMinW;
                        info = 1;
                    }

                    if (xj==Zero) {
                        continue;
                    }

                    if (Tjj<One) {
                        if (xj>bigNum*Tjj) {
                            const ElementType rec = One / xj;
                            x     *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                    }
                    ElementType sr, si;
                    ladiv(x(j1), x(n+j1), tmp, z, sr, si);
                    x(j1)   = sr;
                    x(n+j1) = si;
                    xj = abs(x(j1)) + abs(x(n+j1));
//
//                  Scale x if necessary to avoid overflow when adding a
//                  multiple of column j1 of T.
//
                    if (xj>One) {
                        const ElementType rec = One / xj;
                        if (work(j1)>(bigNum-xMax)*rec) {
                            x     *= rec;
                            scale *= rec;
                        }
                    }

                    if (j1>1) {
                        x(_(1,j1-1))     -= x(j1) * T(_(1,j1-1),j1);
                        x(_(n+1,n+j1-1)) -= x(n+j1) * T(_(1,j1-1),j1);

                        x(1)   += b(j1) * x(n+j1);
                        x(n+1) -= b(j1) * x(j1);

                        xMax = Zero;
                        for (IndexType k=1; k<=j1-1; ++k) {
                            xMax = max(xMax, abs(x(k)) + abs(x(k+n)));
                        }
                    }

                } else {
//
//                  Meet 2 by 2 diagonal block
//
                    D(1,1) = x(j1);
                    D(2,1) = x(j2);
                    D(1,2) = x(n+j1);
                    D(2,2) = x(n+j2);

                    ElementType scaleC;
                    IndexType iErr = laln2(false, IndexType(2), sMinW, One,
                                           T(_(j1,j1+1),_(j1,j1+1)), One, One,
                                           D, Zero, -w, V, scaleC, xNorm);
                    if (iErr!=0) {
                        info = 2;
                    }

                    if (scaleC!=One) {
                        x     *= scaleC;
                        scale *= scaleC;
                    }
                    x(j1) = V(1,1);
                    x(j2) = V(2,1);
                    x(n+j1) = V(1,2);
                    x(n+j2) = V(2,2);
//
//                  Scale X(J1), .... to avoid overflow in
//                  updating right hand side.
//
                    const ElementType xj = max(abs(V(1,1))+abs(V(1,2)),
                                               abs(V(2,1))+abs(V(2,2)));
                    if (xj>One) {
                        const ElementType rec = One / xj;
                        if (max(work(j1), work(j2))>(bigNum-xMax)*rec) {
                            x     *= rec;
                            scale *= rec;
                        }
                    }
//
//                  Update the right-hand side.
//
                    if (j1>1) {
                        x(_(1,j1-1)) -= x(j1) * T(_(1,j1-1),j1);
                        x(_(1,j1-1)) -= x(j2) * T(_(1,j1-1),j2);

                        x(_(n+1,n+j1-1)) -= x(n+j1) * T(_(1,j1-1),j1);
                        x(_(n+1,n+j1-1)) -= x(n+j2) * T(_(1,j1-1),j2);

                        x(1)   += b(j1)*x(n+j1);
                        x(1)   += b(j2)*x(n+j2);
                        x(n+1) -= b(j1)*x(j1);
                        x(n+1) -= b(j2)*x(j2);

                        xMax = Zero;
                        for (IndexType k=1; k<=j1-1; ++k) {
                            xMax = max(abs(x(k)) + abs(x(k+n)), xMax);
                        }
                    }

                }
            }

        } else {
//
//          Solve (T + iB)**T*(p+iq) = c+id
//
            IndexType jNext = 1;
            for (IndexType j=1; j<=n; ++j) {
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
//                  1 by 1 diagonal block
//
//                  Scale if necessary to avoid overflow in forming the
//                  right-hand side element by inner product.
//
                    ElementType xj = abs(x(j1)) + abs(x(j1+n));
                    if (xMax>One) {
                        const ElementType rec = One / xMax;
                        if (work(j1)>(bigNum-xj)*rec) {
                            x     *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                    }

                    x(j1)   -= T(_(1,j1-1),j1) * x(_(1,j1-1));
                    x(n+j1) -= T(_(1,j1-1),j1) * x(_(n+1,n+j1-1));
                    if (j1>1) {
                        x(j1)   -= b(j1) * x(n+1);
                        x(n+j1) += b(j1) * x(1);
                    }
                    xj = abs(x(j1)) + abs(x(j1+n));

                    ElementType z = w;
                    if (j1==1) {
                        z = b(1);
                    }
//
//                  Scale if necessary to avoid overflow in
//                  complex division
//
                    ElementType Tjj = abs(T(j1,j1)) + abs(z);
                    ElementType tmp = T(j1,j1);
                    if (Tjj<sMinW) {
                        tmp = sMinW;
                        Tjj = sMinW;
                        info = 1;
                    }
                    if (Tjj<One) {
                        if (xj>bigNum*Tjj) {
                            const ElementType rec = One / xj;
                            x     *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                    }
                    ElementType sr, si;
                    ladiv(x(j1), x(n+j1), tmp, -z, sr, si);
                    x(j1)   = sr;
                    x(j1+n) = si;
                    xMax = max(abs(x(j1)) + abs(x(j1+n)), xMax);

                } else {
//
//                  2 by 2 diagonal block
//
//                  Scale if necessary to avoid overflow in forming the
//                  right-hand side element by inner product.
//
                    const ElementType xj = max(abs(x(j1)) + abs(x(n+j1)),
                                               abs(x(j2)) + abs(x(n+j2)));
                    if (xMax>One) {
                        const ElementType rec = One / xMax;
                        if (max(work(j1), work(j2))>(bigNum-xj)/xj) {
                            x     *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                    }

                    D(1,1) = x(j1) - T(_(1,j1-1),j1) * x(_(1,j1-1));
                    D(2,1) = x(j2) - T(_(1,j1-1),j2) * x(_(1,j1-1));
                    D(1,2) = x(n+j1) - T(_(1,j1-1),j1) * x(_(n+1,n+j1-1));
                    D(2,2) = x(n+j2) - T(_(1,j1-1),j2) * x(_(n+1,n+j1-1));
                    D(1,1) -= b(j1) * x(n+1);
                    D(2,1) -= b(j2) * x(n+1);
                    D(1,2) += b(j1) * x(1);
                    D(2,2) += b(j2) * x(1);

                    ElementType scaleC;
                    IndexType iErr = laln2(true, IndexType(2), sMinW, One,
                                           T(_(j1,j1+1),_(j1,j1+1)), One, One,
                                           D, Zero, w, V, scaleC, xNorm);
                    if (iErr!=0) {
                        info = 2;
                    }

                    if (scaleC!=One) {
                        x     *= scaleC;
                        scale *= scaleC;
                    }
                    x(j1) = V(1,1);
                    x(j2) = V(2,1);
                    x(n+j1) = V(1,2);
                    x(n+j2) = V(2,2);
                    xMax = max(abs(x(j1))+abs(x(n+j1)),
                               abs(x(j2))+abs(x(n+j2)),
                               xMax);

                }

            }

        }

    }

    return info;
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MT, typename VB, typename W, typename SCALE, typename VX,
          typename VWORK>
typename DenseVector<VX>::IndexType
laqtr_native(bool                  trans,
             bool                  real,
             const GeMatrix<MT>    &T,
             const DenseVector<VB> &b,
             const W               &w,
             SCALE                 &scale,
             DenseVector<VX>       &x,
             DenseVector<VWORK>    &work)
{
    typedef typename DenseVector<VX>::ElementType  ElementType;

    const LOGICAL    LTRAN  = trans;
    const LOGICAL    LREAL  = real;
    const INTEGER    N      = T.numRows();
    const INTEGER    LDT    = T.leadingDimension();
    const DOUBLE     _W     = w;
    DOUBLE           _SCALE = scale;
    INTEGER          INFO;

    if (IsSame<ElementType, DOUBLE>::value) {
        LAPACK_IMPL(dlaqtr)(&LTRAN,
                            &LREAL,
                            &N,
                            T.data(),
                            &LDT,
                            b.data(),
                            &_W,
                            &_SCALE,
                            x.data(),
                            work.data(),
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO>=0);

    scale   = _SCALE;

    return INFO;
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================
template <typename MT, typename VB, typename W, typename SCALE, typename VX,
          typename VWORK>
typename DenseVector<VX>::IndexType
laqtr(bool                  trans,
      bool                  real,
      const GeMatrix<MT>    &T,
      const DenseVector<VB> &b,
      const W               &w,
      SCALE                 &scale,
      DenseVector<VX>       &x,
      DenseVector<VWORK>    &work)
{
    typedef typename DenseVector<VX>::IndexType  IndexType;

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());

    const IndexType n = T.numRows();

    if (!real) {
        ASSERT(b.firstIndex()==1);
        ASSERT(b.length()==n);
    }

    ASSERT(x.length()==2*n);
    ASSERT(work.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    SCALE                                   scale_org  = scale;
    typename DenseVector<VX>::NoView        x_org      = x;
    typename DenseVector<VWORK>::NoView     work_org   = work;
#   endif

    IndexType info = laqtr_generic(trans, real, T, b, w, scale, x, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    SCALE                                   scale_generic  = scale;
    typename DenseVector<VX>::NoView        x_generic      = x;
    typename DenseVector<VWORK>::NoView     work_generic   = work;
//
//  restore output arguments
//
    scale = scale_org;
    x     = x_org;
    work  = work_org;
//
//  Compare results
//
    IndexType _info = laqtr_native(trans, real, T, b, w, scale, x, work);

    bool failed = false;
    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }

    if (! isIdentical(x_generic, x, "x_generic", "x")) {
        std::cerr << "CXXLAPACK: x_generic = " << x_generic << std::endl;
        std::cerr << "F77LAPACK: x = " << x << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info= " << info<< std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MT, typename VB, typename W, typename SCALE, typename VX,
          typename VWORK>
typename VX::IndexType
laqtr(bool          trans,
      bool          real,
      const MT      &T,
      const VB      &b,
      const W       &w,
      SCALE         &scale,
      VX            &x,
      VWORK         &work)
{
    typedef typename VX::IndexType  IndexType;

    CHECKPOINT_ENTER;
    const IndexType info = laqtr(trans, real, T, b, w, scale, x, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LAQTR_TCC
