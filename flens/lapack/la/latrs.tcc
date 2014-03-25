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

/*
       SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
      $                   CNORM, INFO )
       SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
      $                   CNORM, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LATRS_TCC
#define FLENS_LAPACK_LA_LATRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- latrs [real variant] ------------------------------------------------------

template <typename MA, typename VX, typename SCALE, typename CNORM>
typename RestrictTo<IsReal<typename MA::ElementType>::value,
         void>::Type
latrs_impl(Transpose             trans,
           bool                  normIn,
           const TrMatrix<MA>    &A,
           DenseVector<VX>       &x,
           SCALE                 &scale,
           DenseVector<CNORM>    &cNorm)
{
    using std::abs;
    using std::max;
    using std::min;

    typedef typename DenseVector<VX>::ElementType  T;
    typedef typename DenseVector<VX>::IndexType    IndexType;

    const T Zero(0), Half(0.5), One(1);
    const Underscore<IndexType>  _;
    const IndexType n = A.dim();
    const bool upper = (A.upLo()==Upper);
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
//
//  Determine machine dependent parameters to control overflow.
//
    const T smallNum = lamch<T>(SafeMin) / lamch<T>(Precision);
    const T bigNum = One / smallNum;
    scale = One;

    if (!normIn) {
//
//      Compute the 1-norm of each column, not including the diagonal.
//
        if (upper) {
//
//          A is upper triangular.
//
            for (IndexType j=1; j<=n; ++j) {
                cNorm(j) = blas::asum(A(_(1,j-1),j));
            }
        } else {
//
//          A is lower triangular.
//
            for (IndexType j=1; j<=n-1; ++j) {
                cNorm(j) = blas::asum(A(_(j+1,n),j));
            }
            cNorm(n) = Zero;
        }
    }
//
//  Scale the column norms by TSCAL if the maximum element in CNORM is
//  greater than BIGNUM.
//
    const IndexType iMax = blas::iamax(cNorm);
    const T tMax = cNorm(iMax);

    T tScale;
    if (tMax<=bigNum) {
        tScale = One;
    } else {
        tScale = One / (smallNum*tMax);
        cNorm *= tScale;
    }
//
//  Compute a bound on the computed solution vector to see if the
//  Level 2 BLAS routine DTRSV can be used.
//
    IndexType jFirst, jLast, jEnd, jInc;

    const IndexType j = blas::iamax(x);
    T xMax = abs(x(j));

    T xBound = xMax;
    T grow;

    if (trans==NoTrans) {
//
//      Compute the growth in A * x = b.
//
        if (upper) {
            jFirst = n;
            jLast  = 1;
            jInc   = -1;
        } else {
            jFirst = 1;
            jLast  = n;
            jInc   = 1;
        }
        jEnd = jLast + jInc;

        if (tScale!=One) {
            grow = Zero;
        } else {
            if (A.diag()==NonUnit) {
//
//              A is non-unit triangular.
//
//              Compute GROW = 1/G(j) and XBND = 1/M(j).
//              Initially, G(0) = max{x(i), i=1,...,n}.
//
                grow = One / max(xBound, smallNum);
                xBound = grow;
                bool grothFactorTooSmall = false;
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        grothFactorTooSmall = true;
                        break;
                    }
//
//                  M(j) = G(j-1) / abs(A(j,j))
//
                    const T  Tjj = abs(A(j,j));
                    xBound = min(xBound, min(One,Tjj)*grow);
                    if (Tjj+cNorm(j)>=smallNum) {
//
//                      G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
//
                        grow *= Tjj / (Tjj+cNorm(j));
                    } else {
//
//                      G(j) could overflow, set GROW to 0.
//
                        grow = Zero;
                    }
                }
                if (!grothFactorTooSmall) {
                    grow = xBound;
                }
            } else {
//
//              A is unit triangular.
//
//              Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                grow = min(One, One/max(xBound, smallNum));
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        break;
                    }
//
//                  G(j) = G(j-1)*( 1 + CNORM(j) )
//
                    grow *= One / (One+cNorm(j));
                }
            }
        }
    } else {
//
//      Compute the growth in A**T * x = b.
//
        if (upper) {
            jFirst = 1;
            jLast  = n;
            jInc   = 1;
        } else {
            jFirst = n;
            jLast  = 1;
            jInc   = -1;
        }
        jEnd = jLast + jInc;

        if (tScale!=One) {
            grow = Zero;
        } else {
            if (A.diag()==NonUnit) {
//
//              A is non-unit triangular.
//
//              Compute GROW = 1/G(j) and XBND = 1/M(j).
//              Initially, M(0) = max{x(i), i=1,...,n}.
//
                grow = One / max(xBound, smallNum);
                xBound = grow;

                bool grothFactorTooSmall = false;
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        grothFactorTooSmall = false;
                        break;
                    }
//
//                  G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
//
                    const T xj = One + cNorm(j);
                    grow = min(grow, xBound / xj);
//
//                  M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
//
                    const T Tjj = abs(A(j,j));
                    if (xj>Tjj) {
                        xBound *= Tjj/xj;
                    }
                }
                if (!grothFactorTooSmall) {
                    grow = min(grow, xBound);
                }
            } else {
//
//              A is unit triangular.
//
//              Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                grow = min(One, One / max(xBound, smallNum));
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        break;
                    }
//
//                  G(j) = ( 1 + CNORM(j) )*G(j-1)
//
                    const T xj = One + cNorm(j);
                    grow /= xj;
                }
            }
        }
    }

    if ((grow*tScale)>smallNum) {
//
//      Use the Level 2 BLAS solve if the reciprocal of the bound on
//      elements of X is not too small.
//
        blas::sv(trans, A, x);
    } else {
//
//      Use a Level 1 BLAS solve, scaling intermediate results.
//
        if (xMax>bigNum) {
//
//          Scale X so that its components are less than or equal to
//          BIGNUM in absolute value.
//
            scale = bigNum /xMax;
            x *= scale;
            xMax = bigNum;
        }

        if (trans==NoTrans) {
//
//          Solve A * x = b
//
            for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
//
                T xj = abs(x(j));
                T TjjS;

                bool skip = false;

                if (A.diag()==NonUnit) {
                    TjjS = A(j,j) * tScale;
                } else {
                    TjjS = tScale;
                    if (tScale==One) {
                        skip = true;
                    }
                }
                if (!skip) {
                    const T  Tjj = abs(TjjS);
                    if (Tjj>smallNum) {
//
//                      abs(A(j,j)) > SMLNUM:
//
                        if (Tjj<One) {
                            if (xj>Tjj*bigNum) {
//
//                              Scale x by 1/b(j).
//
                                const T rec = One / xj;
                                x *= rec;
                                scale *= rec;
                                xMax *= rec;
                            }
                        }
                        x(j) /= TjjS;
                        xj = abs(x(j));
                    } else if (Tjj>Zero) {
//
//                      0 < abs(A(j,j)) <= SMLNUM:
//
                        if (xj>Tjj*bigNum) {
//
//                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
//                          to avoid overflow when dividing by A(j,j).
//
                            T rec = (Tjj*bigNum) / xj;
                            if (cNorm(j)>One) {
//
//                              Scale by 1/CNORM(j) to avoid overflow when
//                              multiplying x(j) times column j.
//
                                rec /= cNorm(j);
                            }
                            x *= rec;
                            scale *= rec;
                            xMax *= rec;
                        }
                        x(j) /= TjjS;
                        xj = abs(x(j));
                    } else {
//
//                      A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//                      scale = 0, and compute a solution to A*x = 0.
//
                        x = Zero;
                        x(j) = One;
                        xj = One;
                        scale = Zero;
                        xMax = Zero;
                    }
                }
//
//              Scale x if necessary to avoid overflow when adding a
//              multiple of column j of A.
//
                if (xj>One) {
                    T rec = One / xj;
                    if (cNorm(j)>(bigNum-xMax)*rec) {
//
//                      Scale x by 1/(2*abs(x(j))).
//
                        rec *= Half;
                        x *= rec;
                        scale *= rec;
                    }
                } else if (xj*cNorm(j)>(bigNum-xMax)) {
//
//                  Scale x by 1/2.
//
                    x *= Half;
                    scale *= Half;
                }

                if (upper) {
                    if (j>1) {
//
//                      Compute the update
//                      x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
//
                        x(_(1,j-1)) -= (x(j)*tScale) * A(_(1,j-1),j);
                        const IndexType i = blas::iamax(x(_(1,j-1)));
                        xMax = abs(x(i));
                    }
                } else {
                    if (j<n) {
//
//                      Compute the update
//                      x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
//
                        x(_(j+1,n)) -= (x(j)*tScale) * A(_(j+1,n),j);
                        const IndexType i = j + blas::iamax(x(_(j+1,n)));
                        xMax = abs(x(i));
                    }
                }
            }

        } else {
//
//          Solve A**T * x = b
//
            for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//              Compute x(j) = b(j) - sum A(k,j)*x(k).
//                                    k<>j
//
                T xj = abs(x(j));
                T uScale = tScale;
                T rec = One / max(xMax, One);
                T TjjS = Zero;

                if (cNorm(j)>(bigNum-xj)*rec) {
//
//                  If x(j) could overflow, scale x by 1/(2*XMAX).
//
                    rec *= Half;
                    if (A.diag()==NonUnit) {
                        TjjS = A(j,j)*tScale;
                    } else {
                        TjjS = tScale;
                    }
                    const T  Tjj = abs(TjjS);
                    if (Tjj>One) {
//
//                      Divide by A(j,j) when scaling x if A(j,j) > 1.
//
                        rec = min(One, rec*Tjj);
                        uScale /= TjjS;
                    }
                    if (rec<One) {
                        x *= rec;
                        scale *= rec;
                        xMax *= rec;
                    }
                }

                T sumJ = Zero;
                if (uScale==One) {
//
//                  If the scaling needed for A in the dot product is 1,
//                  call DDOT to perform the dot product.
//
                    if (upper) {
                        sumJ = A(_(1,j-1),j) * x(_(1,j-1));
                    } else if (j<n) {
                        sumJ = A(_(j+1,n),j) * x(_(j+1,n));
                    }
                } else {
//
//                  Otherwise, use in-line code for the dot product.
//
                    if (upper) {
                        for (IndexType i=1; i<=j-1; ++i) {
                            sumJ += (A(i,j)*uScale)*x(i);
                        }
                    } else if (j<n) {
                        for (IndexType i=j+1; i<=n; ++i) {
                            sumJ += (A(i,j)*uScale)*x(i);
                        }
                    }
                }
                if (uScale==tScale) {
//
//                  Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
//                  was not used to scale the dotproduct.
//
                    x(j) -= sumJ;
                    xj = abs(x(j));

                    bool skip = false;
                    if (A.diag()==NonUnit) {
                        TjjS = A(j,j)*tScale;
                    } else {
                        TjjS = tScale;
                        if (tScale==One) {
                            skip = true;
                        }
                    }
                    if (!skip) {
//
//                      Compute x(j) = x(j) / A(j,j), scaling if necessary.
//
                        const T  Tjj = abs(TjjS);
                        if (Tjj>smallNum) {
//
//                          abs(A(j,j)) > SMLNUM:
//
                            if (Tjj<One) {
                                if (xj>Tjj*bigNum) {
//
//                                  Scale X by 1/abs(x(j)).
//
                                    rec = One / xj;
                                    x *= rec;
                                    scale *= rec;
                                    xMax *= rec;
                                }
                            }
                            x(j) /= TjjS;
                        } else if (Tjj>Zero) {
//
//                          0 < abs(A(j,j)) <= SMLNUM:
//
                            if (xj>Tjj*bigNum) {
//
//                              Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
//
                                const T rec = (Tjj*bigNum) / xj;
                                x *= rec;
                                scale *= rec;
                                xMax *= rec;
                            }
                            x(j) /= TjjS;
                        } else {
//
//                          A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//                          scale = 0, and compute a solution to A**T*x = 0.
//
                            x = Zero;
                            x(j) = One;
                            scale = Zero;
                            xMax = Zero;
                        }
                    }
                } else {
//
//                  Compute x(j) := x(j) / A(j,j)  - sumj if the dot
//                  product has already been divided by 1/A(j,j).
//
                    x(j) = x(j)/TjjS - sumJ;
                }
                xMax = max(xMax, abs(x(j)));
            }
        }
        scale /= tScale;
    }
//
//  Scale the column norms by 1/TSCAL for return.
//
    if (tScale!=One) {
        cNorm *= One/tScale;
    }
}

//
//  Auxiliary function for complex variant
//
template <typename T>
T
abs2(const std::complex<T> &x)
{
    using std::abs;
    using std::imag;
    using std::real;

    const T Two(2);

    return abs(real(x)/Two) + abs(imag(x)/Two);
}

//
//  Complex variant
//
template <typename MA, typename VX, typename SCALE, typename CNORM>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         void>::Type
latrs_impl(Transpose             trans,
           bool                  normIn,
           const TrMatrix<MA>    &A,
           DenseVector<VX>       &x,
           SCALE                 &scale,
           DenseVector<CNORM>    &cNorm)
{
    using std::abs;
    using std::max;
    using std::min;
    using cxxblas::abs1;

    typedef typename DenseVector<VX>::ElementType    T;
    typedef typename ComplexTrait<T>::PrimitiveType  PT;
    typedef typename DenseVector<VX>::IndexType      IndexType;

    const PT Zero(0), Half(0.5), One(1), Two(2);
    const Underscore<IndexType>  _;
    const IndexType n = A.dim();
    const bool upper = (A.upLo()==Upper);
//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
//
//  Determine machine dependent parameters to control overflow.
//
    PT smallNum = lamch<PT>(SafeMin);
    PT bigNum   = One / smallNum;
    labad(smallNum, bigNum);
    smallNum   /= lamch<PT>(Precision);
    bigNum      = One/smallNum;

    scale = One;

    if (!normIn) {
//
//      Compute the 1-norm of each column, not including the diagonal.
//
        if (upper) {
//
//          A is upper triangular.
//
            for (IndexType j=1; j<=n; ++j) {
                cNorm(j) = blas::asum(A(_(1,j-1),j));
            }
        } else {
//
//          A is lower triangular.
//
            for (IndexType j=1; j<=n-1; ++j) {
                cNorm(j) = blas::asum(A(_(j+1,n),j));
            }
            cNorm(n) = Zero;
        }
    }
//
//  Scale the column norms by TSCAL if the maximum element in CNORM is
//  greater than BIGNUM.
//
    const IndexType iMax = blas::iamax(cNorm);
    const PT        tMax = cNorm(iMax);

    PT tScale;
    if (tMax<=bigNum*Half) {
        tScale = One;
    } else {
        tScale = Half / (smallNum*tMax);
        cNorm *= tScale;
    }
//
//  Compute a bound on the computed solution vector to see if the
//  Level 2 BLAS routine DTRSV can be used.
//
    IndexType jFirst, jLast, jEnd, jInc;

    PT xMax = Zero;
    for (IndexType j=1; j<=n; ++j) {
        xMax = max(xMax, abs2(x(j)));
    }

    PT xBound = xMax;
    PT grow;

    if (trans==NoTrans) {
//
//      Compute the growth in A * x = b.
//
        if (upper) {
            jFirst = n;
            jLast  = 1;
            jInc   = -1;
        } else {
            jFirst = 1;
            jLast  = n;
            jInc   = 1;
        }
        jEnd = jLast + jInc;

        if (tScale!=One) {
            grow = Zero;
        } else {
            if (A.diag()==NonUnit) {
//
//              A is non-unit triangular.
//
//              Compute GROW = 1/G(j) and XBND = 1/M(j).
//              Initially, G(0) = max{x(i), i=1,...,n}.
//
                grow = Half / max(xBound, smallNum);
                xBound = grow;
                bool grothFactorTooSmall = false;
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        grothFactorTooSmall = true;
                        break;
                    }

                    const T  TjjS = A(j,j);
                    const PT Tjj  = abs1(TjjS);

                    if (Tjj>=smallNum) {
//
//                      M(j) = G(j-1) / abs(A(j,j))
                        xBound = min(xBound, min(One, Tjj)*grow);
                    } else {
//
//                      M(j) could overflow, set XBND to 0.
//
                        xBound = Zero;
                    }

                    if (Tjj+cNorm(j)>=smallNum) {
//
//                      G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
//
                        grow *= Tjj / (Tjj+cNorm(j));
                    } else {
//
//                      G(j) could overflow, set GROW to 0.
//
                        grow = Zero;
                    }
                }
                if (!grothFactorTooSmall) {
                    grow = xBound;
                }
            } else {
//
//              A is unit triangular.
//
//              Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                grow = min(One, Half/max(xBound, smallNum));
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        break;
                    }
//
//                  G(j) = G(j-1)*( 1 + CNORM(j) )
//
                    grow *= One / (One+cNorm(j));
                }
            }
        }
    } else {
//
//      Compute the growth in A**T * x = b.
//
        if (upper) {
            jFirst = 1;
            jLast  = n;
            jInc   = 1;
        } else {
            jFirst = n;
            jLast  = 1;
            jInc   = -1;
        }
        jEnd = jLast + jInc;

        if (tScale!=One) {
            grow = Zero;
        } else {
            if (A.diag()==NonUnit) {
//
//              A is non-unit triangular.
//
//              Compute GROW = 1/G(j) and XBND = 1/M(j).
//              Initially, M(0) = max{x(i), i=1,...,n}.
//
                grow = Half / max(xBound, smallNum);
                xBound = grow;

                bool grothFactorTooSmall = false;
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        grothFactorTooSmall = false;
                        break;
                    }
//
//                  G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
//
                    const PT xj = One + cNorm(j);
                    grow = min(grow, xBound / xj);

                    const T  TjjS = A(j,j);
                    const PT Tjj  = abs1(TjjS);

                    if (Tjj>=smallNum) {
//
//                      M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
//
                        if (xj>Tjj) {
                            xBound *= Tjj/xj;
                        }
                    } else {
//
//                      M(j) could overflow, set XBND to 0.
//
                        xBound = Zero;
                    }
                }
                if (!grothFactorTooSmall) {
                    grow = min(grow, xBound);
                }
            } else {
//
//              A is unit triangular.
//
//              Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                grow = min(One, Half / max(xBound, smallNum));
                for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//                  Exit the loop if the growth factor is too small.
//
                    if (grow<=smallNum) {
                        break;
                    }
//
//                  G(j) = ( 1 + CNORM(j) )*G(j-1)
//
                    const PT xj = One + cNorm(j);
                    grow /= xj;
                }
            }
        }
    }

    if ((grow*tScale)>smallNum) {
//
//      Use the Level 2 BLAS solve if the reciprocal of the bound on
//      elements of X is not too small.
//
        blas::sv(trans, A, x);
    } else {
//
//      Use a Level 1 BLAS solve, scaling intermediate results.
//
        if (xMax>bigNum*Half) {
//
//          Scale X so that its components are less than or equal to
//          BIGNUM in absolute value.
//
            scale = (bigNum*Half) /xMax;
            x *= scale;
            xMax = bigNum;
        } else {
            xMax *= Two;
        }

        if (trans==NoTrans) {
//
//          Solve A * x = b
//
            for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
//
                PT xj = abs1(x(j));
                T  TjjS;

                bool skip = false;

                if (A.diag()==NonUnit) {
                    TjjS = A(j,j) * tScale;
                } else {
                    TjjS = tScale;
                    if (tScale==One) {
                        skip = true;
                    }
                }
                if (!skip) {
                    const PT  Tjj = abs1(TjjS);
                    if (Tjj>smallNum) {
//
//                      abs(A(j,j)) > SMLNUM:
//
                        if (Tjj<One) {
                            if (xj>Tjj*bigNum) {
//
//                              Scale x by 1/b(j).
//
                                const PT rec = One / xj;
                                x     *= rec;
                                scale *= rec;
                                xMax  *= rec;
                            }
                        }
                        x(j) = ladiv(x(j), TjjS);
                        xj    = abs1(x(j));
                    } else if (Tjj>Zero) {
//
//                      0 < abs(A(j,j)) <= SMLNUM:
//
                        if (xj>Tjj*bigNum) {
//
//                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
//                          to avoid overflow when dividing by A(j,j).
//
                            PT rec = (Tjj*bigNum) / xj;
                            if (cNorm(j)>One) {
//
//                              Scale by 1/CNORM(j) to avoid overflow when
//                              multiplying x(j) times column j.
//
                                rec /= cNorm(j);
                            }
                            x     *= rec;
                            scale *= rec;
                            xMax  *= rec;
                        }
                        x(j) = ladiv(x(j), TjjS);
                        xj = abs1(x(j));
                    } else {
//
//                      A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//                      scale = 0, and compute a solution to A*x = 0.
//
                        x = Zero;
                        x(j) = One;
                        xj = One;
                        scale = Zero;
                        xMax = Zero;
                    }
                }
//
//              Scale x if necessary to avoid overflow when adding a
//              multiple of column j of A.
//
                if (xj>One) {
                    PT rec = One / xj;
                    if (cNorm(j)>(bigNum-xMax)*rec) {
//
//                      Scale x by 1/(2*abs(x(j))).
//
                        rec   *= Half;
                        x     *= rec;
                        scale *= rec;
                    }
                } else if (xj*cNorm(j)>(bigNum-xMax)) {
//
//                  Scale x by 1/2.
//
                    x     *= Half;
                    scale *= Half;
                }

                if (upper) {
                    if (j>1) {
//
//                      Compute the update
//                      x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
//
                        x(_(1,j-1)) -= (x(j)*tScale) * A(_(1,j-1),j);
                        const IndexType i = blas::iamax(x(_(1,j-1)));
                        xMax = abs1(x(i));
                    }
                } else {
                    if (j<n) {
//
//                      Compute the update
//                      x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
//
                        x(_(j+1,n)) -= (x(j)*tScale) * A(_(j+1,n),j);
                        const IndexType i = j + blas::iamax(x(_(j+1,n)));
                        xMax = abs1(x(i));
                    }
                }
            }

        } else if (trans==Trans) {
//
//          Solve A**T * x = b
//
            for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//              Compute x(j) = b(j) - sum A(k,j)*x(k).
//                                    k<>j
//
                PT xj     = abs1(x(j));
                T  uScale = tScale;
                PT rec    = One / max(xMax, One);
                T  TjjS   = Zero;

                if (cNorm(j)>(bigNum-xj)*rec) {
//
//                  If x(j) could overflow, scale x by 1/(2*XMAX).
//
                    rec *= Half;
                    if (A.diag()==NonUnit) {
                        TjjS = A(j,j)*tScale;
                    } else {
                        TjjS = tScale;
                    }
                    const PT  Tjj = abs1(TjjS);
                    if (Tjj>One) {
//
//                      Divide by A(j,j) when scaling x if A(j,j) > 1.
//
                        rec = min(One, rec*Tjj);
                        uScale = ladiv(uScale, TjjS);
                    }
                    if (rec<One) {
                        x     *= rec;
                        scale *= rec;
                        xMax  *= rec;
                    }
                }

                T sumJ = Zero;
                if (uScale==One) {
//
//                  If the scaling needed for A in the dot product is 1,
//                  call DDOT to perform the dot product.
//
                    if (upper) {
                        sumJ = A(_(1,j-1),j) * x(_(1,j-1));
                    } else if (j<n) {
                        sumJ = A(_(j+1,n),j) * x(_(j+1,n));
                    }
                } else {
//
//                  Otherwise, use in-line code for the dot product.
//
                    if (upper) {
                        for (IndexType i=1; i<=j-1; ++i) {
                            sumJ += (A(i,j)*uScale)*x(i);
                        }
                    } else if (j<n) {
                        for (IndexType i=j+1; i<=n; ++i) {
                            sumJ += (A(i,j)*uScale)*x(i);
                        }
                    }
                }

                if (uScale==T(tScale)) {
//
//                  Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
//                  was not used to scale the dotproduct.
//
                    x(j) -= sumJ;
                    xj = abs1(x(j));

                    bool skip = false;
                    if (A.diag()==NonUnit) {
                        TjjS = A(j,j)*tScale;
                    } else {
                        TjjS = tScale;
                        if (tScale==One) {
                            skip = true;
                        }
                    }
                    if (!skip) {
//
//                      Compute x(j) = x(j) / A(j,j), scaling if necessary.
//
                        const PT  Tjj = abs1(TjjS);
                        if (Tjj>smallNum) {
//
//                          abs(A(j,j)) > SMLNUM:
//
                            if (Tjj<One) {
                                if (xj>Tjj*bigNum) {
//
//                                  Scale X by 1/abs(x(j)).
//
                                    rec = One / xj;
                                    x     *= rec;
                                    scale *= rec;
                                    xMax  *= rec;
                                }
                            }
                            x(j) = ladiv(x(j), TjjS);
                        } else if (Tjj>Zero) {
//
//                          0 < abs(A(j,j)) <= SMLNUM:
//
                            if (xj>Tjj*bigNum) {
//
//                              Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
//
                                const PT rec = (Tjj*bigNum) / xj;
                                x     *= rec;
                                scale *= rec;
                                xMax  *= rec;
                            }
                            x(j) = ladiv(x(j), TjjS);
                        } else {
//
//                          A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//                          scale = 0, and compute a solution to A**T*x = 0.
//
                            x = Zero;
                            x(j) = One;
                            scale = Zero;
                            xMax = Zero;
                        }
                    }
                } else {
//
//                  Compute x(j) := x(j) / A(j,j)  - sumj if the dot
//                  product has already been divided by 1/A(j,j).
//
                    x(j) = ladiv(x(j), TjjS) - sumJ;
                }
                xMax = max(xMax, abs1(x(j)));
            }
        } else {
//
//          Solve A**H * x = b
//
            for (IndexType j=jFirst; j!=jEnd; j+=jInc) {
//
//              Compute x(j) = b(j) - sum A(k,j)*x(k).
//                                   k<>j
//
                PT xj     = abs1(x(j));
                T  uScale = tScale;
                PT rec    = One / max(xMax, One);
                T  TjjS   = Zero;

                if (cNorm(j)>(bigNum-xj)*rec) {
//
//                  If x(j) could overflow, scale x by 1/(2*XMAX).
//
                    rec *= Half;
                    if (A.diag()==NonUnit) {
                        TjjS = conj(A(j,j))*tScale;
                    } else {
                        TjjS = tScale;
                    }
                    const PT  Tjj = abs1(TjjS);
                    if (Tjj>One) {
//
//                      Divide by A(j,j) when scaling x if A(j,j) > 1.
//
                        rec = min(One, rec*Tjj);
                        uScale = ladiv(uScale, TjjS);
                    }
                    if (rec<One) {
                        x     *= rec;
                        scale *= rec;
                        xMax  *= rec;
                    }
                }

                T sumJ = Zero;
                if (uScale==One) {
//
//                  If the scaling needed for A in the dot product is 1,
//                  call ZDOTC to perform the dot product.
//
                    if (upper) {
                        sumJ = conjugate(A(_(1,j-1),j)) * x(_(1,j-1));
                    } else if (j<n) {
                        sumJ = conjugate(A(_(j+1,n),j)) * x(_(j+1,n));
                    }
                } else {
//
//                  Otherwise, use in-line code for the dot product.
//
                    if (upper) {
                        for (IndexType i=1; i<=j-1; ++i) {
                            sumJ += (conj(A(i,j))*uScale)*x(i);
                        }
                    } else if (j<n) {
                        for (IndexType i=j+1; i<=n; ++i) {
                            sumJ += (conj(A(i,j))*uScale)*x(i);
                        }
                    }
                }

                if (uScale==T(tScale)) {
//
//                  Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
//                  was not used to scale the dotproduct.
//
                    x(j) -= sumJ;
                    xj = abs1(x(j));

                    bool skip = false;
                    if (A.diag()==NonUnit) {
                        TjjS = conj(A(j,j))*tScale;
                    } else {
                        TjjS = tScale;
                        if (tScale==One) {
                            skip = true;
                        }
                    }
                    if (!skip) {
//
//                      Compute x(j) = x(j) / A(j,j), scaling if necessary.
//
                        const PT  Tjj = abs1(TjjS);
                        if (Tjj>smallNum) {
//
//                          abs(A(j,j)) > SMLNUM:
//
                            if (Tjj<One) {
                                if (xj>Tjj*bigNum) {
//
//                                  Scale X by 1/abs(x(j)).
//
                                    rec = One / xj;
                                    x     *= rec;
                                    scale *= rec;
                                    xMax  *= rec;
                                }
                            }
                            x(j) = ladiv(x(j), TjjS);
                        } else if (Tjj>Zero) {
//
//                          0 < abs(A(j,j)) <= SMLNUM:
//
                            if (xj>Tjj*bigNum) {
//
//                              Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
//
                                const PT rec = (Tjj*bigNum) / xj;
                                x     *= rec;
                                scale *= rec;
                                xMax  *= rec;
                            }
                            x(j) = ladiv(x(j), TjjS);
                        } else {
//
//                          A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//                          scale = 0 and compute a solution to A**H *x = 0.
//
                            x = Zero;
                            x(j) = One;
                            scale = Zero;
                            xMax = Zero;
                        }
                    } 
                    
                } else {
//
//                  Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
//                  product has already been divided by 1/A(j,j).
//
                    x(j) = ladiv(x(j), TjjS) - sumJ;
                }
                xMax = max(xMax, abs1(x(j)));
                
            }
        }
        scale /= tScale;
    }
//
//  Scale the column norms by 1/TSCAL for return.
//
    if (tScale!=One) {
        cNorm *= One/tScale;
    }
}


} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- latrs [real and complex variant] ------------------------------------------

template <typename MA, typename VX, typename SCALE, typename CNORM>
void
latrs_impl(Transpose             trans,
           bool                  normIn,
           const TrMatrix<MA>    &A,
           DenseVector<VX>       &x,
           SCALE                 &scale,
           DenseVector<CNORM>    &cNorm)
{
    typedef typename TrMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::latrs<IndexType>(getF77Char(A.upLo()),
                                                 getF77Char(trans),
                                                 getF77Char(A.diag()),
                                                 (normIn) ? 'Y' : 'N',
                                                 A.dim(),
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 x.data(),
                                                 scale,
                                                 cNorm.data());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- latrs [real variant] ------------------------------------------------------

template <typename MA, typename VX, typename SCALE, typename VCNORM>
typename RestrictTo<IsRealTrMatrix<MA>::value
                 && IsRealDenseVector<VX>::value
                 && IsNotComplex<SCALE>::value
                 && IsRealDenseVector<VCNORM>::value,
         void>::Type
latrs(Transpose trans,
      bool      normIn,
      const MA  &A,
      VX        &&x,
      SCALE     &scale,
      VCNORM    &&cNorm)
{
    LAPACK_DEBUG_OUT("latrs [real]");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VX>::Type        VectorX;
    typedef typename RemoveRef<SCALE>::Type     Scale;
    typedef typename RemoveRef<VCNORM>::Type    VectorCNorm;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(x.firstIndex()==1);
    ASSERT(x.length()==A.dim());
    ASSERT(cNorm.firstIndex()==1);
    ASSERT(cNorm.length()==A.dim());
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename VectorX::NoView        x_org     = x;
    Scale                           scale_org = scale;
    typename VectorCNorm::NoView    cNorm_org = cNorm;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::latrs_impl(trans, normIn, A, x, scale, cNorm);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorX::NoView     x_generic     = x;
    Scale                        scale_generic = scale;
    typename VectorCNorm::NoView cNorm_generic = cNorm;

    x     = x_org;
    scale = scale_org;
    cNorm = cNorm_org;

    external::latrs_impl(trans, normIn, A, x, scale, cNorm);

    bool failed = false;
    if (! isIdentical(x_generic, x, "x_generic", "x")) {
        std::cerr << "CXXLAPACK: x_generic = " << x_generic << std::endl;
        std::cerr << "F77LAPACK: x = " << x << std::endl;
        failed = true;
    }

    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }

    if (! isIdentical(cNorm_generic, cNorm, "cNorm_generic", "cNorm")) {
        std::cerr << "CXXLAPACK: cNorm_generic = "
                  << cNorm_generic << std::endl;
        std::cerr << "F77LAPACK: cNorm = " << cNorm << std::endl;
        failed = true;
    }


    if (failed) {
        std::cerr << "x_org = " << x_org << std::endl;
        std::cerr << "scale_org = " << scale_org << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: latrs" << std::endl;
    }
#   endif
}

//-- latrs [complex variant] ---------------------------------------------------

template <typename MA, typename VX, typename SCALE, typename VCNORM>
typename RestrictTo<IsComplexTrMatrix<MA>::value
                 && IsComplexDenseVector<VX>::value
                 && IsNotComplex<SCALE>::value
                 && IsRealDenseVector<VCNORM>::value,
         void>::Type
latrs(Transpose trans,
      bool      normIn,
      const MA  &A,
      VX        &&x,
      SCALE     &scale,
      VCNORM    &&cNorm)
{
    LAPACK_DEBUG_OUT("latrs [complex]");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VX>::Type        VectorX;
    typedef typename RemoveRef<SCALE>::Type     Scale;
    typedef typename RemoveRef<VCNORM>::Type    VectorCNorm;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(x.firstIndex()==1);
    ASSERT(x.length()==A.dim());
    ASSERT(cNorm.firstIndex()==1);
    ASSERT(cNorm.length()==A.dim());
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename VectorX::NoView        x_org     = x;
    Scale                           scale_org = scale;
    typename VectorCNorm::NoView    cNorm_org = cNorm;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::latrs_impl(trans, normIn, A, x, scale, cNorm);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorX::NoView     x_generic     = x;
    Scale                        scale_generic = scale;
    typename VectorCNorm::NoView cNorm_generic = cNorm;

    x     = x_org;
    scale = scale_org;
    cNorm = cNorm_org;

    external::latrs_impl(trans, normIn, A, x, scale, cNorm);

    bool failed = false;
    if (! isIdentical(x_generic, x, "x_generic", "x")) {
        std::cerr << "CXXLAPACK: x_generic = " << x_generic << std::endl;
        std::cerr << "F77LAPACK: x = " << x << std::endl;
        failed = true;
    }

    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }

    if (! isIdentical(cNorm_generic, cNorm, "cNorm_generic", "cNorm")) {
        std::cerr << "CXXLAPACK: cNorm_generic = "
                  << cNorm_generic << std::endl;
        std::cerr << "F77LAPACK: cNorm = " << cNorm << std::endl;
        failed = true;
    }


    if (failed) {
        if(trans==NoTrans) {
            std::cerr << "trans = NoTrans" << std::endl; 
        } else if (trans==Conj) {
            std::cerr << "trans = Conj" << std::endl; 
        } else if (trans==Trans) {
            std::cerr << "trans = Trans" << std::endl; 
        } else if (trans==ConjTrans) {
            std::cerr << "trans = ConjTrans" << std::endl; 
        }
        std::cerr << "normIn = " << normIn << std::endl;
        std::cerr << "x_org = " << x_org << std::endl;
        std::cerr << "scale_org = " << scale_org << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: latrs" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LATRS_TCC