/*
 *   Copyright (c) 2011, Iris Haecker
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
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *   
 *   $COPYRIGHT$
 *   
 *   Additional copyrights may follow
 *   
 *   $HEADER$
 *   
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *   
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *     
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *     
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
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

#ifndef CXXLAPACK_AUX_LASET_TCC
#define CXXLAPACK_AUX_LASET_TCC 1

#include <cxxblas/cxxblas.h>
#include <cxxlapack/cxxlapack.h>

namespace cxxlapack {

template <typename IndexType, typename T>
IndexType
latrs(StorageOrder order,
      StorageUpLo upLo, Transpose trans, Diag diag, bool normIn,
      IndexType n, const T *A, IndexType ldA, T *x, T &scale, T *cNorm)
{
    assert(order==ColMajor);

    IndexType info = 0;
//
//  Test the input parameters.
//
    if (n<0) {
        info = -5;
    } else if (ldA<max(1,n)) {
        info = -7;
    }
    if (info!=0) {
        // xerbla
        return info;
    }
//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Determine machine dependent parameters to control overflow.
//
    const T smallNum = lamch<T>(SafeMin) / lamch<T>(Precision);
    const T bigNum = T(1) / smallNum;
    scale = 1;

    if (!normIn) {
//
//      Compute the 1-norm of each column, not including the diagonal.
//
        if (upLo==Upper) {
//
//          A is upper triangular.
//
            cNorm[0] = 0;
            for (IndexType j=1; j<n; ++j) {
                cNorm[j] = asum(j, A+j*ldA, 1)
            }
        } else {
//
//          A is lower triangular.
//
            for (IndexType j=0; j<n-1; ++j) {
                cNorm[j] = asum(n-1-j, A+(j+1)*ldA+1, 1);
            }
            cNorm[n-1] = 0;
        }
//
//      Scale the column norms by tScale if the maximum element in cNorm is
//      greater than bigNum.
//
        IndexType iMax = iamax(n, cNorm, IndexType(1));
        T         tMax = cnorm[iMax];
        T         tScale;

        if (tMax<=bigNum) {
            tScale = 1;
        } else {
            tScale = T(1) / (smallNum*tMax);
            scal(n, tScale, cNorm, 1);
        }
//
//      Compute a bound on the computed solution vector to see if the
//      Level 2 BLAS routine STRSV can be used.
//
        IndexType j = iamax(n, x, 1);
        T         xMax = abs(x[j]);
        T         xBound = xMax;

        if (trans==NoTrans) {
//
//          Compute the growth in A * x = b.
//
            if (upLo==Upper) {
                jFirst = n-1;
                jLast = 0;
                jInc = -1;
            } else {
                jFirst = 0;
                jLast = n-1
                jInc = 1;
            }

            if (tScale!=T(1)) {
                grow = 0;
            } else {
                if (diag==NonUnit) {
//
//                  A is non-unit triangular.
//
//                  Compute GROW = 1/G(j) and XBND = 1/M(j).
//                  Initially, G(0) = max{x(i), i=1,...,n}.
//
                    grow = T(1) / max(xBound, smallNum);
                    xBound = grow;
                    for (IndexType j=jFirst; j!=jLast; j+=jInc) {
//
//                      Exit the loop if the growth factor is too small.
//
                        if (grow<=smallNum) {
                            break;
                        }
//
//                      M(j) = G(j-1) / abs(A(j,j))
//
                        tjj = abs(A[j+ldA*j]);
                        xBound = min(xBound, min(T(1), tjj)*grow);
                        if (tjj+cNorm[j]>=smallNum) {
//
//                          G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
//
                            grow = grow*(tjj/(tjj+cNorm[j]));
                        } else {
//
//                          G(j) could overflow, set GROW to 0.
//
                            grow = 0;
                        }
                    }
                    grow = xBound;
                } else {
//
//                  A is unit triangular.
//
//                  Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                    grow = min(T(1), T(1) / max(xBound, smallNum));
                    for (IndexType j=jFirst; j!=jLast; j+=jInc) {
//
//                      Exit the loop if the growth factor is too small.
//
                        if (grow<=smallNum) {
                            break;
                        }
//
//                      G(j) = G(j-1)*( 1 + CNORM(j) )
//
                        grow = grow*(T(1) / (T(1)+cNorm[j]));
                    }
                }
            }
        } else {
//
//          Compute the growth in A**T * x = b.
//
            if (upLo==Upper) {
                jFirst = 0;
                jLast = n-1;
                jInc = 1;
            } else {
                jFirst = n-1;
                jLast = 0;
                jInc = -1;
            }

            if (tScale!=T(1)) {
                grow = T(0);
            } else {
                if (diag==NonUnit) {
//
//                  A is non-unit triangular.
//
//                  Compute GROW = 1/G(j) and XBND = 1/M(j).
//                  Initially, M(0) = max{x(i), i=1,...,n}.
//
                    grow = T(1) / max(xBound, smallNum);
                    xBound = grow;

                    for (IndexType j=jFirst; j!=jLast; j+=jInc) {
//
//                      Exit the loop if the growth factor is too small.
//
                        if (grow<=smallNum) {
                            break;
                        }
//
//                      G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
//
                        xj = T(1) + cNorm[j];
                        grow = min(grow, xBound / xj);
//
//                      M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
//
                        tjj = abs(A[j*(ldA+1)]);
                        if (xj>tjj) {
                            xBound = xBound*(tjj/xj);
                        }
                    }
                    grow = min(grow, xBound);
                } else {
//
//                  A is unit triangular.
//
//                  Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
//
                    grow = min(T(1), T(1) / max(xBound, smallNum));
                    for (IndexType j=jFirst; j!=jLast; j+=jInc) {
//
//                      Exit the loop if the growth factor is too small.
//
                        if (grow<=smallNum) {
                            break;
                        }
//
//                      G(j) = ( 1 + CNORM(j) )*G(j-1)
//
                        xj = T(1) + cNorm[j];
                        grow = grow / xj;
                    }
                }
            }
        }
    *
          IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
    *
    *        Use the Level 2 BLAS solve if the reciprocal of the bound on
    *        elements of X is not too small.
    *
             CALL STRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
          ELSE
    *
    *        Use a Level 1 BLAS solve, scaling intermediate results.
    *
             IF( XMAX.GT.BIGNUM ) THEN
    *
    *           Scale X so that its components are less than or equal to
    *           BIGNUM in absolute value.
    *
                SCALE = BIGNUM / XMAX
                CALL SSCAL( N, SCALE, X, 1 )
                XMAX = BIGNUM
             END IF
    *
             IF( NOTRAN ) THEN
    *
    *           Solve A * x = b
    *
                DO 100 J = JFIRST, JLAST, JINC
    *
    *              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
    *
                   XJ = ABS( X( J ) )
                   IF( NOUNIT ) THEN
                      TJJS = A( J, J )*TSCAL
                   ELSE
                      TJJS = TSCAL
                      IF( TSCAL.EQ.ONE )
         $               GO TO 95
                   END IF
                      TJJ = ABS( TJJS )
                      IF( TJJ.GT.SMLNUM ) THEN
    *
    *                    abs(A(j,j)) > SMLNUM:
    *
                         IF( TJJ.LT.ONE ) THEN
                            IF( XJ.GT.TJJ*BIGNUM ) THEN
    *
    *                          Scale x by 1/b(j).
    *
                               REC = ONE / XJ
                               CALL SSCAL( N, REC, X, 1 )
                               SCALE = SCALE*REC
                               XMAX = XMAX*REC
                            END IF
                         END IF
                         X( J ) = X( J ) / TJJS
                         XJ = ABS( X( J ) )
                      ELSE IF( TJJ.GT.ZERO ) THEN
    *
    *                    0 < abs(A(j,j)) <= SMLNUM:
    *
                         IF( XJ.GT.TJJ*BIGNUM ) THEN
    *
    *                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
    *                       to avoid overflow when dividing by A(j,j).
    *
                            REC = ( TJJ*BIGNUM ) / XJ
                            IF( CNORM( J ).GT.ONE ) THEN
    *
    *                          Scale by 1/CNORM(j) to avoid overflow when
    *                          multiplying x(j) times column j.
    *
                               REC = REC / CNORM( J )
                            END IF
                            CALL SSCAL( N, REC, X, 1 )
                            SCALE = SCALE*REC
                            XMAX = XMAX*REC
                         END IF
                         X( J ) = X( J ) / TJJS
                         XJ = ABS( X( J ) )
                      ELSE
    *
    *                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
    *                    scale = 0, and compute a solution to A*x = 0.
    *
                         DO 90 I = 1, N
                            X( I ) = ZERO
       90                CONTINUE
                         X( J ) = ONE
                         XJ = ONE
                         SCALE = ZERO
                         XMAX = ZERO
                      END IF
       95          CONTINUE
    *
    *              Scale x if necessary to avoid overflow when adding a
    *              multiple of column j of A.
    *
                   IF( XJ.GT.ONE ) THEN
                      REC = ONE / XJ
                      IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
    *
    *                    Scale x by 1/(2*abs(x(j))).
    *
                         REC = REC*HALF
                         CALL SSCAL( N, REC, X, 1 )
                         SCALE = SCALE*REC
                      END IF
                   ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
    *
    *                 Scale x by 1/2.
    *
                      CALL SSCAL( N, HALF, X, 1 )
                      SCALE = SCALE*HALF
                   END IF
    *
                   IF( UPPER ) THEN
                      IF( J.GT.1 ) THEN
    *
    *                    Compute the update
    *                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
    *
                         CALL SAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X,
         $                           1 )
                         I = ISAMAX( J-1, X, 1 )
                         XMAX = ABS( X( I ) )
                      END IF
                   ELSE
                      IF( J.LT.N ) THEN
    *
    *                    Compute the update
    *                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
    *
                         CALL SAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1,
         $                           X( J+1 ), 1 )
                         I = J + ISAMAX( N-J, X( J+1 ), 1 )
                         XMAX = ABS( X( I ) )
                      END IF
                   END IF
      100       CONTINUE
    *
             ELSE
    *
    *           Solve A**T * x = b
    *
                DO 140 J = JFIRST, JLAST, JINC
    *
    *              Compute x(j) = b(j) - sum A(k,j)*x(k).
    *                                    k<>j
    *
                   XJ = ABS( X( J ) )
                   USCAL = TSCAL
                   REC = ONE / MAX( XMAX, ONE )
                   IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
    *
    *                 If x(j) could overflow, scale x by 1/(2*XMAX).
    *
                      REC = REC*HALF
                      IF( NOUNIT ) THEN
                         TJJS = A( J, J )*TSCAL
                      ELSE
                         TJJS = TSCAL
                      END IF
                         TJJ = ABS( TJJS )
                         IF( TJJ.GT.ONE ) THEN
    *
    *                       Divide by A(j,j) when scaling x if A(j,j) > 1.
    *
                            REC = MIN( ONE, REC*TJJ )
                            USCAL = USCAL / TJJS
                         END IF
                      IF( REC.LT.ONE ) THEN
                         CALL SSCAL( N, REC, X, 1 )
                         SCALE = SCALE*REC
                         XMAX = XMAX*REC
                      END IF
                   END IF
    *
                   SUMJ = ZERO
                   IF( USCAL.EQ.ONE ) THEN
    *
    *                 If the scaling needed for A in the dot product is 1,
    *                 call SDOT to perform the dot product.
    *
                      IF( UPPER ) THEN
                         SUMJ = SDOT( J-1, A( 1, J ), 1, X, 1 )
                      ELSE IF( J.LT.N ) THEN
                         SUMJ = SDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                      END IF
                   ELSE
    *
    *                 Otherwise, use in-line code for the dot product.
    *
                      IF( UPPER ) THEN
                         DO 110 I = 1, J - 1
                            SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
      110                CONTINUE
                      ELSE IF( J.LT.N ) THEN
                         DO 120 I = J + 1, N
                            SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
      120                CONTINUE
                      END IF
                   END IF
    *
                   IF( USCAL.EQ.TSCAL ) THEN
    *
    *                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
    *                 was not used to scale the dotproduct.
    *
                      X( J ) = X( J ) - SUMJ
                      XJ = ABS( X( J ) )
                      IF( NOUNIT ) THEN
                         TJJS = A( J, J )*TSCAL
                      ELSE
                         TJJS = TSCAL
                         IF( TSCAL.EQ.ONE )
         $                  GO TO 135
                      END IF
    *
    *                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
    *
                         TJJ = ABS( TJJS )
                         IF( TJJ.GT.SMLNUM ) THEN
    *
    *                       abs(A(j,j)) > SMLNUM:
    *
                            IF( TJJ.LT.ONE ) THEN
                               IF( XJ.GT.TJJ*BIGNUM ) THEN
    *
    *                             Scale X by 1/abs(x(j)).
    *
                                  REC = ONE / XJ
                                  CALL SSCAL( N, REC, X, 1 )
                                  SCALE = SCALE*REC
                                  XMAX = XMAX*REC
                               END IF
                            END IF
                            X( J ) = X( J ) / TJJS
                         ELSE IF( TJJ.GT.ZERO ) THEN
    *
    *                       0 < abs(A(j,j)) <= SMLNUM:
    *
                            IF( XJ.GT.TJJ*BIGNUM ) THEN
    *
    *                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
    *
                               REC = ( TJJ*BIGNUM ) / XJ
                               CALL SSCAL( N, REC, X, 1 )
                               SCALE = SCALE*REC
                               XMAX = XMAX*REC
                            END IF
                            X( J ) = X( J ) / TJJS
                         ELSE
    *
    *                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
    *                       scale = 0, and compute a solution to A**T*x = 0.
    *
                            DO 130 I = 1, N
                               X( I ) = ZERO
      130                   CONTINUE
                            X( J ) = ONE
                            SCALE = ZERO
                            XMAX = ZERO
                         END IF
      135             CONTINUE
                   ELSE
    *
    *                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
    *                 product has already been divided by 1/A(j,j).
    *
                      X( J ) = X( J ) / TJJS - SUMJ
                   END IF
                   XMAX = MAX( XMAX, ABS( X( J ) ) )
      140       CONTINUE
             END IF
             SCALE = SCALE / TSCAL
          END IF
    *
    *     Scale the column norms by 1/TSCAL for return.
    *
          IF( TSCAL.NE.ONE ) THEN
             CALL SSCAL( N, ONE / TSCAL, CNORM, 1 )
          END IF
    *
          RETURN
    *
    *     End of SLATRS
    *
          END


} // namespace cxxlapack

#endif // CXXLAPACK_AUX_LASET_TCC
