/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,
      $                   SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                  --
 *
 *  -- Contributed by Zlatko Drmac of the University of Zagreb and     --
 *  -- Kresimir Veselic of the Fernuniversitaet Hagen                  --
 *  -- April 2011                                                      --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 * This routine is also part of SIGMA (version 1.23, October 23. 2008.)
 * SIGMA is a library of algorithms for highly accurate algorithms for
 * computation of SVD, PSVD, QSVD, (H,K)-SVD, and for solution of the
 * eigenvalue problems Hx = lambda M x, H M x = lambda x with H, M > 0.
 *
 */

#ifndef FLENS_LAPACK_SVD_SVJ0_TCC
#define FLENS_LAPACK_SVD_SVJ0_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj0_generic(SVJ::JobV                                  jobV,
             GeMatrix<MA>                               &A,
             DenseVector<VD>                            &d,
             DenseVector<VSVA>                          &sva,
             GeMatrix<MV>                               &V,
             const typename GeMatrix<MA>::ElementType   &eps,
             const typename GeMatrix<MA>::ElementType   &safeMin,
             const typename GeMatrix<MA>::ElementType   &tol,
             typename GeMatrix<MA>::IndexType           nSweep,
             DenseVector<VWORK>                         &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType  Zero(0), Half(0.5), One(1);

    const Underscore<IndexType>  _;

    ElementType fastr_data[5];
    DenseVectorView<ElementType>
        fastr  = typename DenseVectorView<ElementType>::Engine(5, fastr_data);

    const IndexType  m     = A.numRows();
    const IndexType  n     = A.numCols();

    auto _work = work(_(1,m));

    const bool applyV   = (jobV==SVJ::ApplyV);
    const bool rhsVec   = (jobV==SVJ::ComputeV) || applyV;

    const ElementType rootEps = sqrt(eps);
    const ElementType rootSafeMin = sqrt(safeMin);
    const ElementType small = safeMin / eps;
    const ElementType big = One / safeMin;
    const ElementType rootBig = One / rootSafeMin;
    const ElementType bigTheta = One / rootEps;
    const ElementType rootTol = sqrt(tol);

    IndexType  info = 0;
//
//  -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#-
//
    const IndexType emptsw = n*(n-1)/2;
    IndexType notRot = 0;
    fastr(1) = Zero;
//
//  -#- Row-cyclic pivot strategy with de Rijk's pivoting -#-
//

    IndexType swBand = 0;
//[TP] SWBAND is a tuning parameter. It is meaningful and effective
//  if SGESVJ is used as a computational routine in the preconditioned
//  Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
//  ......

    const IndexType kbl = min(IndexType(8),n);
//[TP] KBL is a tuning parameter that defines the tile size in the
//  tiling of the p-q loops of pivot pairs. In general, an optimal
//  value of KBL depends on the matrix dimensions and on the
//  parameters of the computer's memory.
//
    IndexType nbl = n / kbl;
    if (nbl*kbl!=n) {
        ++nbl;
    }

    const IndexType blSkip = pow(kbl,2) + 1;
//[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

    const IndexType rowSkip = min(IndexType(5),kbl);
//[TP] ROWSKIP is a tuning parameter.

    const IndexType lkAhead = 1;
//[TP] LKAHEAD is a tuning parameter.
    IndexType pSkipped = 0;

    ElementType aapp, aapp0, aaqq, aapq, aqoap, apoaq;
    ElementType tmp, cs, sn, t, theta, thetaSign;

    bool rotOk;
    bool converged = false;

    for (IndexType i=1; i<=nSweep; ++i) {
//  .. go go go ...
//
       ElementType max_aapq = Zero;
       ElementType max_sinj = Zero;

       IndexType iswRot = 0;

       pSkipped = 0;
       notRot = 0;

       for (IndexType ibr=1; ibr<=nbl; ++ibr) {

          IndexType igl = (ibr-1)*kbl + 1;

          for (IndexType ir1=0; ir1<=min(lkAhead, nbl-ibr); ++ir1) {

             igl += ir1*kbl;

             for (IndexType p=igl; p<=min(igl+kbl-1,n-1); ++p) {

//  .. de Rijk's pivoting
                IndexType q = blas::iamax(sva(_(p,n))) + p - 1;
                if (p!=q) {
                   blas::swap(A(_,p), A(_,q));
                   if (rhsVec) {
                      blas::swap(V(_,p), V(_,q));
                   }
                   swap(sva(p),sva(q));
                   swap(d(p),d(q));
                }

                if (ir1==0) {
//
//     Column norms are periodically updated by explicit
//     norm computation.
//     Caveat:
//     Some BLAS implementations compute DNRM2(M,A(1,p),1)
//     as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may result in
//     overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and
//     undeflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
//     Hence, DNRM2 cannot be trusted, not even in the case when
//     the true norm is far from the under(over)flow boundaries.
//     If properly implemented DNRM2 is available, the IF-THEN-ELSE
//     below should read "AAPP = DNRM2( M, A(1,p), 1 ) * D(p)".
//
                   if (sva(p)<rootBig && sva(p)>rootSafeMin) {
                      sva(p) = blas::nrm2(A(_,p))*d(p);
                   } else {
                      ElementType tmp = Zero;
                      aapp = One;
                      lassq(A(_,p), tmp, aapp);
                      sva(p) = tmp*sqrt(aapp)*d(p);
                   }
                   aapp = sva(p);
                } else {
                   aapp = sva(p);
                }

                if (aapp>Zero) {

                   pSkipped = 0;

                   for (IndexType q=p+1; q<=min(igl+kbl-1,n); ++q) {

                      aaqq = sva(q);

                      if (aaqq>Zero) {

                         aapp0 = aapp;
                         if (aaqq>=One) {
                            rotOk = small*aapp<=aaqq;
                            if (aapp<big/aaqq) {
                               aapq = A(_,p)*A(_,q)*d(p)*d(q)/aaqq/aapp;
                            } else {
                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, d(p), _work);
                               aapq = _work*A(_,q)*d(q)/aaqq;
                            }
                         } else {
                            rotOk = aapp<=aaqq/small;
                            if (aapp>small/aaqq) {
                               aapq = A(_,p)*A(_,q)*d(p)*d(q)/aaqq/aapp;
                            } else {
                               _work = A(_,q);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, d(q), _work);
                               aapq = _work*A(_,p)*d(p)/aapp;
                            }
                         }
//
                         max_aapq = max(max_aapq, abs(aapq));
//
//     TO rotate or NOT to rotate, THAT is the question ...
//
                         if (abs(aapq)>tol) {
//
//        .. rotate
//        ROTATED = ROTATED + ONE
//
                            if (ir1==0) {
                               notRot = 0;
                               pSkipped = 0;
                               ++iswRot;
                            }

                            if (rotOk) {

                               aqoap = aaqq / aapp;
                               apoaq = aapp / aaqq;
                               theta = -Half*abs(aqoap-apoaq)/aapq;

                               if (abs(theta)>bigTheta) {

                                  t = Half / theta;
                                  fastr(3) =  t*d(p)/d(q);
                                  fastr(4) = -t*d(q)/d(p);
                                  blas::rotm(A(_,p), A(_,q), fastr);
                                  if (rhsVec) {
                                     blas::rotm(V(_,p), V(_,q), fastr);
                                  }
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));
                                  max_sinj = max(max_sinj, abs(t));

                               } else {
//
//              .. choose correct signum for THETA and rotate
//
                                  thetaSign = -sign(One, aapq);
                                  t = One
                                     / (theta+thetaSign*sqrt(One+theta*theta));
                                  cs = sqrt(One/(One+t*t));
                                  sn = t*cs;

                                  max_sinj = max(max_sinj, abs(sn));
                                  sva(q) = aaqq
                                          * sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));

                                  apoaq = d(p) / d(q);
                                  aqoap = d(q) / d(p);
                                  if (d(p)>=One) {
                                     if (d(q)>=One) {
                                        fastr(3) =  t*apoaq;
                                        fastr(4) = -t*aqoap;
                                        d(p) *= cs;
                                        d(q) *= cs;
                                        blas::rotm(A(_,p), A(_,q), fastr);
                                        if (rhsVec) {
                                           blas::rotm(V(_,p), V(_,q), fastr);
                                        }
                                     } else {
                                        A(_,p) -= t*aqoap*A(_,q);
                                        A(_,q) += cs*sn*apoaq*A(_,p);
                                        d(p) *= cs;
                                        d(q) /= cs;
                                        if (rhsVec) {
                                           V(_,p) -= t*aqoap*V(_,q);
                                           V(_,q) += cs*sn*apoaq*V(_,p);
                                        }
                                     }
                                  } else {
                                     if (d(q)>=One) {
                                        A(_,q) += t*apoaq*A(_,p);
                                        A(_,p) -= cs*sn*aqoap*A(_,q);
                                        d(p) /= cs;
                                        d(q) *= cs;
                                        if (rhsVec) {
                                           V(_,q) += t*apoaq*V(_,p);
                                           V(_,p) -= cs*sn*aqoap*V(_,q);
                                        }
                                     } else {
                                        if (d(p)>=d(q)) {
                                           A(_,p) -= t*aqoap*A(_,q);
                                           A(_,q) += cs*sn*apoaq*A(_,p);
                                           d(p) *= cs;
                                           d(q) /= cs;
                                           if (rhsVec) {
                                              V(_,p) -= t*aqoap*V(_,q);
                                              V(_,q) += cs*sn*apoaq*V(_,p);
                                           }
                                        } else {
                                           A(_,q) += t*apoaq*A(_,p);
                                           A(_,p) -= cs*sn*aqoap*A(_,q);
                                           d(p) /= cs;
                                           d(q) *= cs;
                                           if (rhsVec) {
                                              V(_,q) += t*apoaq*V(_,p);
                                              V(_,p) -= cs*sn*aqoap*V(_,q);
                                           }
                                        }
                                     }
                                  }
                               }
//
                            } else {
//           .. have to use modified Gram-Schmidt like transformation
                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, One, _work);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, One, A(_,q));
                               tmp = -aapq*d(p)/d(q);
                               A(_,q) += tmp*_work;
                               lascl(LASCL::FullMatrix, 0, 0,
                                     One, aaqq, A(_,q));
                               sva(q) = aaqq*sqrt(max(Zero,One-aapq*aapq));
                               max_sinj = max(max_sinj, safeMin);
                            }
//        END IF ROTOK THEN ... ELSE
//
//        In the case of cancellation in updating SVA(q), SVA(p)
//        recompute SVA(q), SVA(p).
                            if (pow(sva(q)/aaqq,2)<=rootEps) {
                               if (aaqq<rootBig && aaqq>rootSafeMin) {
                                  sva(q) = blas::nrm2(A(_,q))*d(q);
                               } else {
                                  t = Zero;
                                  aaqq = One;
                                  lassq(A(_,q), t, aaqq);
                                  sva(q) = t*sqrt(aaqq)*d(q);
                               }
                            }
                            if (aapp/aapp0<=rootEps) {
                               if (aapp<rootBig && aapp>rootSafeMin) {
                                  aapp = blas::nrm2(A(_,p))*d(p);
                               } else {
                                  t = Zero;
                                  aapp = One;
                                  lassq(A(_,p), t, aapp);
                                  aapp = t*sqrt(aapp)*d(p);
                               }
                               sva(p) = aapp;
                            }
//
                         } else {
//     A(:,p) and A(:,q) already numerically orthogonal
                            if (ir1==0) {
                               ++notRot;
                            }
                            ++pSkipped;
                         }
                      } else {
//     A(:,q) is zero column
                         if (ir1==0) {
                            ++notRot;
                         }
                         ++pSkipped;
                      }

                      if (i<=swBand && pSkipped>rowSkip) {
                         if (ir1==0) {
                            aapp = -aapp;
                         }
                         notRot = 0;
                         break;
                      }

                   }
//  END q-LOOP
//

//  bailed out of q-loop

                   sva(p) = aapp;

                } else {
                   sva(p) = aapp;
                   if (ir1==0 && aapp==Zero) {
                      notRot += min(igl+kbl-1,n) - p;
                   }
                }

             }
//  end of the p-loop
//  end of doing the block ( ibr, ibr )
          }
//  end of ir1-loop
//
//........................................................
//... go to the off diagonal blocks
//
          igl = (ibr-1)*kbl + 1;

          for (IndexType jbc=ibr+1; jbc<=nbl; ++jbc) {

             IndexType jgl = (jbc-1)*kbl + 1;
//
//     doing the block at ( ibr, jbc )
//
             IndexType ijblsk = 0;
             for (IndexType p=igl; p<=min(igl+kbl-1,n); ++p) {

                aapp = sva(p);

                if (aapp>Zero) {

                   pSkipped = 0;

                   for (IndexType q=jgl; q<=min(jgl+kbl-1,n); ++q) {

                      aaqq = sva(q);

                      if (aaqq>Zero) {
                         aapp0 = aapp;
//
//  -#- M x 2 Jacobi SVD -#-
//
//     -#- Safe Gram matrix computation -#-
//
                         if (aaqq>=One) {
                            if (aapp>=aaqq) {
                               rotOk = small*aapp<=aaqq;
                            } else {
                               rotOk = small*aaqq<=aapp;
                            }
                            if (aapp<big/aaqq) {
                               aapq = A(_,p)*A(_,q)*d(p)*d(q)/aaqq/aapp;
                            } else {
                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, d(p), _work);
                               aapq = _work*A(_,q)*d(q)/aaqq;
                            }
                         } else {
                            if (aapp>=aaqq) {
                              rotOk = aapp<=aaqq/small;
                            } else {
                              rotOk = aaqq<=aapp/small;
                            }
                            if (aapp>small/aaqq) {
                               aapq = A(_,p)*A(_,q)*d(p)*d(q)/aaqq/aapp;
                            } else {
                               _work = A(_,q);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, d(q), _work);
                               aapq = _work*A(_,p)*d(p)/aapp;
                            }
                         }

                         max_aapq = max(max_aapq, abs(aapq));
//
//     TO rotate or NOT to rotate, THAT is the question ...
//
                         if (abs(aapq)>tol) {
                            notRot = 0;
//        ROTATED  = ROTATED + 1
                            pSkipped = 0;
                            ++iswRot;

                            if (rotOk) {

                               aqoap = aaqq / aapp;
                               apoaq = aapp / aaqq;
                               theta = -Half*abs(aqoap-apoaq)/aapq;
                               if (aaqq>aapp0) {
                                  theta = -theta;
                               }

                               if (abs(theta)>bigTheta) {
                                  t = Half / theta;
                                  fastr(3) =  t*d(p)/d(q);
                                  fastr(4) = -t*d(q)/d(p);
                                  blas::rotm(A(_,p), A(_,q), fastr);
                                  if (rhsVec) {
                                     blas::rotm(V(_,p), V(_,q), fastr);
                                  }
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));
                                  max_sinj = max(max_sinj, abs(t));
                               } else {
//
//              .. choose correct signum for THETA and rotate
//
                                  thetaSign = -sign(One, aapq);
                                  if (aaqq>aapp0) {
                                     thetaSign = -thetaSign;
                                  }
                                  t = One
                                     / (theta+thetaSign*sqrt(One+theta*theta));
                                  cs = sqrt(One/(One+t*t));
                                  sn = t*cs;
                                  max_sinj = max(max_sinj, abs(sn));
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));

                                  apoaq = d(p) / d(q);
                                  aqoap = d(q) / d(p);
                                  if (d(p)>=One) {

                                     if (d(q)>=One) {
                                        fastr(3) =  t*apoaq;
                                        fastr(4) = -t*aqoap;
                                        d(p) *= cs;
                                        d(q) *= cs;
                                        blas::rotm(A(_,p), A(_,q), fastr);
                                        if (rhsVec) {
                                            blas::rotm(V(_,p), V(_,q), fastr);
                                        }
                                     } else {
                                        A(_,p) -= t*aqoap*A(_,q);
                                        A(_,q) += cs*sn*apoaq*A(_,p);
                                        if (rhsVec) {
                                           V(_,p) -= t*aqoap*V(_,q);
                                           V(_,q) += cs*sn*apoaq*V(_,p);
                                        }
                                        d(p) *= cs;
                                        d(q) /= cs;
                                     }
                                  } else {
                                     if (d(q)>=One) {
                                        A(_,q) += t*apoaq*A(_,p);
                                        A(_,p) -= cs*sn*aqoap*A(_,q);
                                        if (rhsVec) {
                                           V(_,q) += t*apoaq*V(_,p);
                                           V(_,p) -= cs*sn*aqoap*V(_,q);
                                        }
                                        d(p) /= cs;
                                        d(q) *= cs;
                                     } else {
                                        if (d(p)>=d(q)) {
                                           A(_,p) -= t*aqoap*A(_,q);
                                           A(_,q) += cs*sn*apoaq*A(_,p);
                                           d(p) *= cs;
                                           d(q) /= cs;
                                           if (rhsVec) {
                                              V(_,p) -= t*aqoap*V(_,q);
                                              V(_,q) += cs*sn*apoaq*V(_,p);
                                           }
                                        } else {
                                           A(_,q) += t*apoaq*A(_,p);
                                           A(_,p) -= cs*sn*aqoap*A(_,q);
                                           d(p) /= cs;
                                           d(q) *= cs;
                                           if (rhsVec) {
                                              V(_,q) += t*apoaq*V(_,p);
                                              V(_,p) -= cs*sn*aqoap*V(_,q);
                                           }
                                        }
                                     }
                                  }
                               }

                            } else {
                               if (aapp>aaqq) {
                                  _work = A(_,p);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aapp, One, _work);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aaqq, One, A(_,q));
                                  tmp = -aapq*d(p)/d(q);
                                  A(_,q) += tmp*_work;
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        One, aaqq, A(_,q));
                                  sva(q) = aaqq*sqrt(max(Zero, One-aapq*aapq));
                                  max_sinj = max(max_sinj, safeMin);
                               } else {
                                  _work = A(_,q);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aaqq, One, _work);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aapp, One, A(_,p));
                                  tmp = -aapq*d(q)/d(p);
                                  A(_,p) += tmp*_work;
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        One, aapp, A(_,p));
                                  sva(p) = aapp*sqrt(max(Zero,One-aapq*aapq));
                                  max_sinj = max(max_sinj, safeMin);
                               }
                            }
//        END IF ROTOK THEN ... ELSE
//
//        In the case of cancellation in updating SVA(q)
//        .. recompute SVA(q)
                            if (pow(sva(q)/aaqq, 2)<=rootEps) {
                               if (aaqq<rootBig && aaqq>rootSafeMin) {
                                  sva(q) = blas::nrm2(A(_,q))*d(q);
                               } else {
                                  t = Zero;
                                  aaqq = One;
                                  lassq(A(_,q), t, aaqq);
                                  sva(q) = t*sqrt(aaqq)*d(q);
                               }
                            }
                            if (pow(aapp/aapp0, 2)<=rootEps) {
                               if (aapp<rootBig && aapp>rootSafeMin) {
                                  aapp = blas::nrm2(A(_,p))*d(p);
                               } else {
                                  t = Zero;
                                  aapp = One;
                                  lassq(A(_,p), t, aapp);
                                  aapp = t*sqrt(aapp)*d(p);
                               }
                               sva(p) = aapp;
                            }
//           end of OK rotation
                         } else {
                            ++notRot;
                            ++pSkipped;
                            ++ijblsk;
                         }
                      } else {
                         ++notRot;
                         ++pSkipped;
                         ++ijblsk;
                      }

                      if (i<=swBand && ijblsk>=blSkip) {
                         sva(p) = aapp;
                         notRot = 0;
                         goto jbcLoopExit;
                      }
                      if (i<=swBand && pSkipped>rowSkip) {
                         aapp = -aapp;
                         notRot = 0;
                         break;
                      }

                   }
//     end of the q-loop

                   sva(p) = aapp;

                } else {
                   if (aapp==Zero) {
                      notRot += min(jgl+kbl-1, n) -jgl + 1;
                   }
                   if (aapp<Zero) {
                      notRot = 0;
                   }
                }

             }
//  end of the p-loop
          }
//  end of the jbc-loop
       jbcLoopExit:
//2011 bailed out of the jbc-loop
          for (IndexType p=igl; p<=min(igl+kbl-1,n); ++p) {
             sva(p) = abs(sva(p));
          }

       }
//2000 :: end of the ibr-loop
//
//  .. update SVA(N)
       if (sva(n)<rootBig  && sva(n)>rootSafeMin) {
          sva(n) = blas::nrm2(A(_,n))*d(n);
       } else {
          ElementType t = Zero;
          aapp = One;
          lassq(A(_,n), t, aapp);
          sva(n) = t*sqrt(aapp)*d(n);
       }
//
//  Additional steering devices
//
       if (i<swBand && (max_aapq<=rootTol || iswRot<=n)) {
           swBand = i;
       }

       if (i>swBand+1 && max_aapq<n*tol && n*max_aapq*max_sinj<tol) {
           converged = true;
           break;
       }

       if (notRot>=emptsw) {
           converged = true;
           break;
       }

    }
//    end i=1:NSWEEP loop
    if (converged) {
//#:) Reaching this point means that during the i-th sweep all pivots were
//    below the given tolerance, causing early exit.
       info = 0;
    } else {
//#:) Reaching this point means that the procedure has comleted the given
//    number of iterations.
        info = nSweep - 1;
    }
//
//  Sort the vector D.
    for (IndexType p=1; p<=n-1; ++p) {
        IndexType q = blas::iamax(sva(_(p,n))) + p - 1;
        if (p!=q) {
            swap(sva(p), sva(q));
            swap(d(p), d(q));
            blas::swap(A(_,p), A(_,q));
            if (rhsVec) {
                blas::swap(V(_,p), V(_,q));
            }
        }
    }
    return info;
}

//== interface for native lapack ===============================================
#ifdef CHECK_CXXLAPACK

template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj0_native(SVJ::JobV                                  jobV,
            GeMatrix<MA>                               &A,
            DenseVector<VD>                            &d,
            DenseVector<VSVA>                          &sva,
            GeMatrix<MV>                               &V,
            const typename GeMatrix<MA>::ElementType   &eps,
            const typename GeMatrix<MA>::ElementType   &safeMin,
            const typename GeMatrix<MA>::ElementType   &tol,
            typename GeMatrix<MA>::IndexType           nSweep,
            DenseVector<VWORK>                         &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const char       JOBV = char(jobV);
    const INTEGER    M = A.numRows();
    const INTEGER    N = A.numCols();
    const INTEGER    LDA = A.leadingDimension();
    const INTEGER    _MV = V.numRows();
    const INTEGER    LDV = V.leadingDimension();
    const DOUBLE     EPS = eps;
    const DOUBLE     SFMIN = safeMin;
    const DOUBLE     TOL = tol;
    const INTEGER    NSWEEP = nSweep;
    const INTEGER    LWORK = work.length();
    INTEGER          INFO;

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_DECL(dgsvj0)(&JOBV,
                            &M,
                            &N,
                            A.data(),
                            &LDA,
                            d.data(),
                            sva.data(),
                            &_MV,
                            V.data(),
                            &LDV,
                            &EPS,
                            &SFMIN,
                            &TOL,
                            &NSWEEP,
                            work.data(),
                            &LWORK,
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO>=0);
    return INFO;
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj0(SVJ::JobV                                  jobV,
     GeMatrix<MA>                               &A,
     DenseVector<VD>                            &d,
     DenseVector<VSVA>                          &sva,
     GeMatrix<MV>                               &V,
     const typename GeMatrix<MA>::ElementType   &eps,
     const typename GeMatrix<MA>::ElementType   &safeMin,
     const typename GeMatrix<MA>::ElementType   &tol,
     typename GeMatrix<MA>::IndexType           nSweep,
     DenseVector<VWORK>                         &work)
{
    using std::max;
    using std::min;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(m>=n);

    ASSERT(d.firstIndex()==1);
    ASSERT(d.length()==n);

    ASSERT(sva.firstIndex()==1);
    ASSERT(sva.length()==n);

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);

    if (jobV==SVJ::ComputeV) {
        ASSERT(V.numCols()==n);
        ASSERT(V.numRows()==n);
    }
    if (jobV==SVJ::ApplyV) {
        ASSERT(V.numCols()==n);
    }

    if (work.length()>0) {
        ASSERT(work.length()>=m);
    }
#   endif
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org     = A;
    typename DenseVector<VD>::NoView    d_org     = d;
    typename DenseVector<VSVA>::NoView  sva_org   = sva;
    typename GeMatrix<MV>::NoView       V_org     = V;
    typename DenseVector<VWORK>::NoView work_org  = work;
#   endif

//
//  Call implementation
//
    IndexType info = svj0_generic(jobV, A, d, sva, V, eps, safeMin, tol,
                                  nSweep, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic     = A;
    typename DenseVector<VD>::NoView    d_generic     = d;
    typename DenseVector<VSVA>::NoView  sva_generic   = sva;
    typename GeMatrix<MV>::NoView       V_generic     = V;
    typename DenseVector<VWORK>::NoView work_generic  = work;
//
//  restore output arguments
//
    A    = A_org;
    d    = d_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = svj0_native(jobV, A, d, sva, V, eps, safeMin, tol,
                                  nSweep, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(d_generic, d, "d_generic", "d")) {
        std::cerr << "CXXLAPACK: d_generic = " << d_generic << std::endl;
        std::cerr << "F77LAPACK: d = " << d << std::endl;
        failed = true;
    }
    if (! isIdentical(sva_generic, sva, "sva_generic", "sva")) {
        std::cerr << "CXXLAPACK: sva_generic = " << sva_generic << std::endl;
        std::cerr << "F77LAPACK: sva = " << sva << std::endl;
        failed = true;
    }
    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
       failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }
    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: svj0.tcc" << std::endl;
        ASSERT(0);
    } else {
        //std::cerr << "passed: svj0.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename MA::IndexType
svj0(SVJ::JobV                                  jobV,
     MA                                         &&A,
     VD                                         &&d,
     VSVA                                       &&sva,
     MV                                         &&V,
     const typename MA::ElementType             &eps,
     const typename MA::ElementType             &safeMin,
     const typename MA::ElementType             &tol,
     typename MA::IndexType                     nSweep,
     VWORK                                      &&work)
{
    typename MA::IndexType   info;

    CHECKPOINT_ENTER;
    info = svj0(jobV, A, d, sva, V, eps, safeMin, tol, nSweep, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SVD_SVJ0_TCC
