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
       SUBROUTINE DGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V,
      $                   LDV, WORK, LWORK, INFO )
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

#ifndef FLENS_LAPACK_GE_SVJ_TCC
#define FLENS_LAPACK_GE_SVJ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj_impl(SVJ::TypeA                typeA,
         SVJ::JobU                 jobU,
         SVJ::JobV                 jobV,
         GeMatrix<MA>              &A,
         DenseVector<VSVA>         &sva,
         GeMatrix<MV>              &V,
         DenseVector<VWORK>        &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using flens::pow;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType  Zero(0), Half(0.5), One(1);
    const IndexType    nSweep = 30;

    const Underscore<IndexType>  _;

    ElementType fastr_data[5];
    DenseVectorView<ElementType>
        fastr  = typename DenseVectorView<ElementType>::Engine(5, fastr_data);

    const IndexType  m     = A.numRows();
    const IndexType  n     = A.numCols();
    const IndexType  mv    = V.numRows();
    const IndexType  lWork = work.length();

    const bool lower    = (typeA==SVJ::Lower);
    const bool upper    = (typeA==SVJ::Upper);
    const bool lhsVec   = (jobU==SVJ::ComputeU);
    const bool controlU = (jobU==SVJ::ControlU);
    const bool applyV   = (jobV==SVJ::ApplyV);
    const bool rhsVec   = (jobV==SVJ::ComputeV) || applyV;

    IndexType  info = 0;
//
//#:) Quick return for void matrix
//
    if ((m==0) || (n==0)) {
        return info;
    }
//
//  Set numerical parameters
//  The stopping criterion for Jacobi rotations is
//
//  max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS
//
//  where EPS is the round-off and CTOL is defined as follows:
//
    ElementType  cTol;
    if (controlU) {
//      ... user controlled
        cTol = work(1);
    } else {
//      ... default
        if (lhsVec || rhsVec) {
            cTol = sqrt(ElementType(m));
        } else {
            cTol = ElementType(m);
        }
    }
//    ... and the machine dependent parameters are
//[!]  (Make sure that DLAMCH() works properly on the target machine.)
//
    const ElementType eps = lamch<ElementType>(Eps);
    const ElementType rootEps = sqrt(eps);
    const ElementType safeMin = lamch<ElementType>(SafeMin);
    const ElementType rootSafeMin = sqrt(safeMin);
    const ElementType small = safeMin / eps;
    const ElementType big = lamch<ElementType>(OverflowThreshold);
//  const ElementType big = One / safeMin;
    const ElementType rootBig = One / rootSafeMin;
    const ElementType bigTheta = One / rootEps;

    const ElementType tol = cTol*eps;
    const ElementType rootTol = sqrt(tol);

    if (ElementType(m)*eps>=One) {
//      we return -4 to keep it compatible with the LAPACK implementation
        return -4;
    }
//
//    Initialize the right singular vector matrix.
//
    if (rhsVec && (!applyV)) {
        V         = Zero;
        V.diag(0) = One;
    }
//
//  Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
//(!) If necessary, scale A to protect the largest singular value
//  from overflow. It is possible that saving the largest singular
//  value destroys the information about the small ones.
//  This initial scaling is almost minimal in the sense that the
//  goal is to make sure that no column norm overflows, and that
//  DSQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
//  in A are detected, the procedure returns with INFO=-6.
//
    ElementType skl = One / sqrt(ElementType(m)*ElementType(n));
    bool noScale = true;
    bool goScale = true;

    ElementType aapp, aaqq, tmp;

    if (lower) {
//      the input matrix is M-by-N lower triangular (trapezoidal)
        for (IndexType p=1; p<=n; ++p) {
            aapp = Zero;
            aaqq = One;
            lassq(A(_(p,m),p), aapp, aaqq);
            if (aapp>big) {
                return -6;
            }
            aaqq = sqrt(aaqq);
            if ((aapp<(big/aaqq)) && noScale) {
                sva(p) = aapp*aaqq;
            } else {
                noScale = false;
                sva(p) = aapp*(aaqq*skl);
                if (goScale) {
                    goScale = false;
                    sva(_(1,p-1)) *= skl;
                }
            }
        }
    } else if (upper) {
//      the input matrix is M-by-N upper triangular (trapezoidal)
        for (IndexType p=1; p<=n; ++p) {
            aapp = Zero;
            aaqq = One;
            lassq(A(_(1,p), p), aapp, aaqq);
            if (aapp>big) {
                return -6;
            }
            aaqq = sqrt(aaqq);
            if ((aapp<(big/aaqq)) && noScale) {
                sva(p) = aapp*aaqq;
            } else {
                noScale = false;
                sva(p) = aapp*(aaqq*skl);
                if (goScale) {
                    goScale = false;
                    sva(_(1,p-1)) *= skl;
                }
            }
        }
    } else {
//      the input matrix is M-by-N general dense
        for (IndexType p=1; p<=n; ++p) {
            aapp = Zero;
            aaqq = One;
            lassq(A(_,p), aapp, aaqq);
            if (aapp>big) {
                return -6;
            }
            aaqq = sqrt(aaqq);
            if ((aapp<(big/aaqq)) && noScale) {
                sva(p) = aapp*aaqq;
            } else {
                noScale = false;
                sva(p) = aapp*(aaqq*skl);
                if (goScale) {
                    goScale = false;
                    sva(_(1,p-1)) *= skl;
                }
            }
        }
    }

    if (noScale) {
        skl = One;
    }
//
//  Move the smaller part of the spectrum from the underflow threshold
//(!) Start by determining the position of the nonzero entries of the
//  array SVA() relative to ( SFMIN, BIG ).
//
    aapp = Zero;
    aaqq = big;
    for (IndexType p=1; p<=n; ++p) {
        if (sva(p)!=Zero) {
            aaqq = min(aaqq,sva(p));
        }
        aapp = max(aapp,sva(p));
    }
//
//#:) Quick return for zero matrix
//
    if (aapp==Zero) {
        if (lhsVec) {
            A = Zero;
            A.diag(0) = One;
        }
        work(1)      = One;
        work(_(2,6)) = Zero;
        return info;
    }
//
//#:) Quick return for one-column matrix
//
    if (n==1) {
        if (lhsVec) {
            lascl(LASCL::FullMatrix, 0, 0, sva(1), skl, A);
        }
        work(1) = One / skl;
        if (sva(1)>=safeMin) {
            work(2) = One;
        } else {
            work(2) = Zero;
        }
        work(_(3,6)) = Zero;
        return info;
    }
//
//  Protect small singular values from underflow, and try to
//  avoid underflows/overflows in computing Jacobi rotations.
//
    ElementType cs, sn;

    sn = sqrt(safeMin/eps);
    tmp = sqrt(big/ElementType(n));
    if ((aapp<=sn) || (aaqq>=tmp) || (sn<=aaqq && aapp<=tmp)) {
        tmp = min(big, tmp/aapp);
//      aaqq  = aaqq*tmp
//      aapp  = aapp*tmp
    } else if ((aaqq<=sn) && (aapp<=tmp)) {
        tmp = min(sn/aaqq, big/(aapp*sqrt(ElementType(n))));
//      aaqq  = aaqq*tmp
//      aapp  = aapp*tmp
    } else if ((aaqq>=sn) && (aapp>=tmp)) {
        tmp = max(sn/aaqq, tmp/aapp);
//      aaqq  = aaqq*tmp
//      aapp  = aapp*tmp
    } else if ((aaqq<=sn) && (aapp>=tmp)) {
        tmp = min(sn/aaqq, big/(sqrt(ElementType(n))*aapp));
//      aaqq  = aaqq*tmp
//      aapp  = aapp*tmp
    } else {
        tmp = One;
    }
//
//  Scale, if necessary
//
    if (tmp!=One) {
        lascl(LASCL::FullMatrix, 0, 0, One, tmp, sva);
    }
    skl *= tmp;
    if (skl!=One) {
        if (upper) {
            lascl(LASCL::UpperTriangular, 0, 0, One, skl, A);
        } else if (lower) {
            lascl(LASCL::LowerTriangular, 0, 0, One, skl, A);
        } else {
            lascl(LASCL::FullMatrix, 0, 0, One, skl, A);
        }
        skl = One / skl;
    }
//
//  Row-cyclic Jacobi SVD algorithm with column pivoting
//
    const IndexType emptsw = n*(n-1)/2;
    IndexType notRot = 0;
    fastr(1) = Zero;
//
//  A is represented in factored form A = A * diag(WORK), where diag(WORK)
//  is initialized to identity. WORK is updated during fast scaled
//  rotations.
//
    work(_(1,n)) = One;
//
//
    IndexType swBand = 3;
//[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
//  if DGESVJ is used as a computational routine in the preconditioned
//  Jacobi SVD algorithm DGESVJ. For sweeps i=1:SWBAND the procedure
//  works on pivots inside a band-like region around the diagonal.
//  The boundaries are determined dynamically, based on the number of
//  pivots above a threshold.
//
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

    const IndexType blSkip = pow(kbl, 2);
//[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
//
    const IndexType rowSkip = min(IndexType(5), kbl);
//[TP] ROWSKIP is a tuning parameter.
//
    const IndexType lkAhead = 1;
//[TP] LKAHEAD is a tuning parameter.
//
//  Quasi block transformations, using the lower (upper) triangular
//  structure of the input matrix. The quasi-block-cycling usually
//  invokes cubic convergence. Big part of this cycle is done inside
//  canonical subspaces of dimensions less than M.
//
    if ((lower || upper) && (n>max(IndexType(64), 4*kbl))) {
//[TP]  The number of partition levels and the actual partition are
//      tuning parameters.
        const IndexType n4 = n / 4;
        const IndexType n2 = n / 2;
        const IndexType n34 = 3*n4;
        const IndexType q = (applyV) ? 0 : 1;
        const IndexType qm = (applyV) ? mv : 0;

        auto _work = work(_(n+1,lWork));

        if (lower) {
//
//          This works very well on lower triangular matrices, in particular
//          in the framework of the preconditioned Jacobi SVD (xGEJSV).
//          The idea is simple:
//          [+ 0 0 0]   Note that Jacobi transformations of [0 0]
//          [+ + 0 0]                                       [0 0]
//          [+ + x 0]   actually work on [x 0]              [x 0]
//          [+ + x x]                    [x x].             [x x]
//
            auto A1 = A(_(    1,m),_(    1, n4));
            auto A2 = A(_( n4+1,m),_( n4+1, n2));
            auto A3 = A(_( n2+1,m),_( n2+1,n34));
            auto A4 = A(_(n34+1,m),_(n34+1,  n));

            auto A12 = A(_(   1,m),_(   1,n2));
            auto A34 = A(_(n2+1,m),_(n2+1, n));

            auto V1 = V(_(      1, n4*q+qm),_(    1,n4));
            auto V2 = V(_( n4*q+1, n2*q+qm),_( n4+1,n2));
            auto V3 = V(_( n2*q+1,n34*q+qm),_( n2+1,n34));
            auto V4 = V(_(n34*q+1, mv*q+qm),_(n34+1,n));

            auto V12 = V(_(      1,n2*q+qm),_(   1,n2));
            auto V34 = V(_( n2*q+1,mv*q+qm),_(n2+1,n));

            auto d1 = work(_(    1, n4));
            auto d2 = work(_( n4+1, n2));
            auto d3 = work(_( n2+1,n34));
            auto d4 = work(_(n34+1,  n));

            auto d12 = work(_(    1,n2));
            auto d34 = work(_( n2+1, n));

            auto sva1 = sva(_(    1, n4));
            auto sva2 = sva(_( n4+1, n2));
            auto sva3 = sva(_( n2+1,n34));
            auto sva4 = sva(_(n34+1,  n));

            auto sva12 = sva(_(   1,n2));
            auto sva34 = sva(_(n2+1, n));

            svj0(jobV, A4, d4, sva4, V4, eps, safeMin, tol, 2, _work);
            svj0(jobV, A3, d3, sva3, V3, eps, safeMin, tol, 2, _work);
            svj1(jobV, n4, A34, d34, sva34, V34, eps, safeMin, tol, 1, _work);
            svj0(jobV, A2, d2, sva2, V2, eps, safeMin, tol, 1, _work);
            svj0(jobV, A1, d1, sva1, V1, eps, safeMin, tol, 1, _work);
            svj1(jobV, n4, A12, d12, sva12, V12, eps, safeMin, tol, 1, _work);

        } else if (upper) {

            auto A1 = A(_(1,   n4),_(    1,   n4));
            auto A2 = A(_(1,   n2),_( n4+1,n4+n4));
            auto A3 = A(_(1,n2+n4),_( n2+1,n2+n4));

            auto A12 = A(_(1,n2),_(1,n2));

            auto V1 = V(_(      1,     n4*q+qm),_(    1,   n4));
            auto V2 = V(_( n4*q+1,(n4+n4)*q+qm),_( n4+1,n4+n4));
            auto V3 = V(_( n2*q+1,(n2+n4)*q+qm),_( n2+1,n2+n4));

            auto V12 = V(_(1,n2*q+qm),_(1,n2));

            auto d1 = work(_(    1,   n4));
            auto d2 = work(_( n4+1,n4+n4));
            auto d3 = work(_( n2+1,n2+n4));

            auto d12 = work(_(1,n2));

            auto sva1 = sva(_(    1,   n4));
            auto sva2 = sva(_( n4+1,n4+n4));
            auto sva3 = sva(_( n2+1,n2+n4));

            auto sva12 = sva(_(1,n2));

            svj0(jobV, A1, d1, sva1, V1, eps, safeMin, tol, 2, _work);
            svj0(jobV, A2, d2, sva2, V2, eps, safeMin, tol, 1, _work);
            svj1(jobV, n4, A12, d12, sva12, V12, eps, safeMin, tol, 1, _work);
            svj0(jobV, A3, d3, sva3, V3, eps, safeMin, tol, 1, _work);

        }

    }
//
//  .. Row-cyclic pivot strategy with de Rijk's pivoting ..
//
    ElementType max_aapq, aapq, aapp0, aqoap, apoaq;
    ElementType max_sinj;
    ElementType theta, thetaSign, t;

    IndexType   i, iswRot, pSkipped, ibr, igl, jgl, ir1, p, q, jbc, ijblsk;

    bool        converged = false;
    bool        rotOk;

    for (i=1; i<=nSweep; ++i) {
//
//  .. go go go ...
//
       max_aapq = Zero;
       max_sinj = Zero;

       iswRot   = 0;
       pSkipped = 0;

       notRot   = 0;
//
//     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
//     1 <= p < q <= N. This is the first step toward a blocked implementation
//     of the rotations. New implementation, based on block transformations,
//     is under development.
//
       for (ibr=1; ibr<=nbl; ++ibr) {

          igl = (ibr-1)*kbl + 1;

          for (ir1=0; ir1<=min(lkAhead,nbl-ibr); ++ir1) {

             igl += ir1*kbl;

             for (p=igl; p<=min(igl+kbl-1,n-1); ++p) {
//
//    .. de Rijk's pivoting
//
                q = blas::iamax(sva(_(p,n))) + p - 1;
                if (p!=q) {
                    blas::swap(A(_,p),A(_,q));
                    if (rhsVec) {
                        blas::swap(V(_,p),V(_,q));
                    }
                    swap(sva(p),sva(q));
                    swap(work(p),work(q));
                }

                if (ir1==0) {
//
//     Column norms are periodically updated by explicit
//     norm computation.
//     Caveat:
//     Unfortunately, some BLAS implementations compute DNRM2(M,A(1,p),1)
//     as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to
//     overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and to
//     underflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
//     Hence, DNRM2 cannot be trusted, not even in the case when
//     the true norm is far from the under(over)flow boundaries.
//     If properly implemented DNRM2 is available, the IF-THEN-ELSE
//     below should read "AAPP = DNRM2( M, A(1,p), 1 ) * WORK(p)".
//

                   if ((sva(p)<rootBig) && (sva(p)>rootSafeMin)) {
                      sva(p) = blas::nrm2(A(_,p)) * work(p);
                   } else {
                       tmp = Zero;
                       aapp = One;
                       lassq(A(_,p), tmp, aapp);
                       sva(p) = tmp * sqrt(aapp) * work(p);
                   }
                   aapp = sva(p);
                } else {
                   aapp = sva(p);
                }

                if (aapp>Zero) {

                   pSkipped = 0;

                   for (q=p+1; q<=min(igl+kbl-1,n); ++q) {

                      aaqq = sva(q);

                      if (aaqq>Zero) {

                         aapp0 = aapp;
                         if (aaqq>=One) {
                            rotOk = (small*aapp<=aaqq);
                            if (aapp<big/aaqq) {
                               aapq = (A(_,p)*A(_,q)*work(p)*work(q)/aaqq)
                                    / aapp;
                            } else {
                               auto _work = work(_(n+1,n+m));

                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, work(p), _work);
                               aapq = _work*A(_,q)*work(q)/aaqq;
                            }
                         } else {
                            rotOk = (aapp<=aaqq/small);
                            if (aapp>(small/aaqq)) {
                                aapq = (A(_,p)*A(_,q)*work(p)*work(q)/aaqq)
                                     / aapp;
                            } else {
                               auto _work = work(_(n+1,n+m));

                               _work = A(_,q);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, work(q), _work);
                               aapq = _work*A(_,p)*work(p)/aapp;
                            }
                         }

                         max_aapq = max(max_aapq, abs(aapq));
//
//     TO rotate or NOT to rotate, THAT is the question ...
//
                         if (abs(aapq)>tol) {
//
//        .. rotate
//[RTD]    ROTATED = ROTATED + ONE
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
                                  fastr(3) =  t*work(p) / work(q);
                                  fastr(4) = -t*work(q) / work(p);
                                  blas::rotm(A(_,p), A(_,q), fastr);
                                  if (rhsVec) {
                                      blas::rotm(V(_,p), V(_,q), fastr);
                                  }
                                  sva(q) =aaqq*sqrt(max(Zero,One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));
                                  max_sinj = max(max_sinj, abs(t));

                               } else {
//
//              .. choose correct signum for THETA and rotate
//
                                  thetaSign = -sign(One,aapq);
                                  t = One
                                    / (theta +thetaSign*sqrt(One+theta*theta));
                                  cs = sqrt(One / (One+t*t));
                                  sn = t*cs;

                                  max_sinj = max(max_sinj, abs(sn));
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));

                                  apoaq = work(p) / work(q);
                                  aqoap = work(q) / work(p);
                                  if (work(p)>=One) {
                                     if (work(q)>=One) {
                                        fastr(3) =  t*apoaq;
                                        fastr(4) = -t*aqoap;
                                        work(p) *= cs;
                                        work(q) *= cs;
                                        blas::rotm(A(_,p), A(_,q), fastr);
                                        if (rhsVec) {
                                            blas::rotm(V(_,p), V(_,q), fastr);
                                        }
                                     } else {
                                        A(_,p) -= t*aqoap*A(_,q);
                                        A(_,q) += cs*sn*apoaq*A(_,p);
                                        work(p) *= cs;
                                        work(q) /= cs;
                                        if (rhsVec) {
                                           V(_,p) -= t*aqoap*V(_,q);
                                           V(_,q) += cs*sn*apoaq*V(_,p);
                                        }
                                     }
                                  } else {
                                     if (work(q)>=One) {
                                        A(_,q) += t*apoaq*A(_,p);
                                        A(_,p) -= cs*sn*aqoap*A(_,q);
                                        work(p) /= cs;
                                        work(q) *= cs;
                                        if (rhsVec) {
                                           V(_,q) += t*apoaq*V(_,p);
                                           V(_,p) -= cs*sn*aqoap*V(_,q);
                                        }
                                     } else {
                                        if (work(p)>=work(q)) {
                                           A(_,p) -= t*aqoap*A(_,q);
                                           A(_,q) += cs*sn*apoaq*A(_,p);
                                           work(p) *= cs;
                                           work(q) /= cs;
                                           if (rhsVec) {
                                              V(_,p) -= t*aqoap*V(_,q);
                                              V(_,q) += cs*sn*apoaq*V(_,p);
                                           }
                                        } else {
                                           A(_,q) += t*apoaq*A(_,p);
                                           A(_,p) -= cs*sn*aqoap*A(_,q);
                                           work(p) /= cs;
                                           work(q) *= cs;
                                           if (rhsVec) {
                                              V(_,q) += t*apoaq*V(_,p);
                                              V(_,p) -= cs*sn*aqoap*V(_,q);
                                           }
                                        }
                                     }
                                  }
                               }

                            } else {
//             .. have to use modified Gram-Schmidt like transformation
                               auto _work = work(_(n+1,lWork));

                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, One, _work);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, One, A(_,q));
                               tmp = -aapq*work(p)/work(q);
                               A(_,q) += tmp*_work;
                               lascl(LASCL::FullMatrix, 0, 0,
                                     One, aaqq, A(_,q));
                               sva(q) = aaqq*sqrt(max(Zero, One-aapq*aapq));
                               max_sinj = max(max_sinj, safeMin);
                            }
//        END IF ROTOK THEN ... ELSE
//
//        In the case of cancellation in updating SVA(q), SVA(p)
//        recompute SVA(q), SVA(p).
//
                            if (pow(sva(q)/aaqq,2)<=rootEps) {
                               if ((aaqq<rootBig) && (aaqq>rootSafeMin)) {
                                  sva(q) = blas::nrm2(A(_,q))*work(q);
                               } else {
                                  t = Zero;
                                  aaqq = One;
                                  lassq(A(_,q), t, aaqq);
                                  sva(q) = t*sqrt(aaqq)*work(q);
                               }
                            }
                            if (aapp/aapp0<=rootEps) {
                               if ((aapp<rootBig) && (aapp>rootSafeMin)) {
                                  aapp = blas::nrm2(A(_,p))*work(p);
                               } else {
                                  t = Zero;
                                  aapp = One;
                                  lassq(A(_,p), t, aapp);
                                  aapp = t*sqrt(aapp)*work(p);
                               }
                               sva(p) = aapp;
                            }
//
                         } else {
//     A(:,p) and A(:,q) already numerically orthogonal
                            if (ir1==0) {
                                ++notRot;
                            }
//[RTD]    SKIPPED  = SKIPPED  + 1
                            ++pSkipped;
                         }
                      } else {
//     A(:,q) is zero column
                         if (ir1==0) {
                             ++notRot;
                         }
                         ++pSkipped;
                      }
//
                      if ((i<=swBand) && (pSkipped>rowSkip)) {
                         if (ir1==0) {
                             aapp = -aapp;
                         }
                         notRot = 0;
                         break;
                      }

                   }
//  END q-LOOP
//
                   sva(p) = aapp;

                } else {
                   sva(p) = aapp;
                   if ((ir1==0) && (aapp==Zero)) {
                       notRot += min(igl+kbl-1,n) - p;
                   }
                }

             }
//  end of the p-loop
//  end of doing the block ( ibr, ibr )
          }
//  end of ir1-loop
//
//... go to the off diagonal blocks
//
          igl = (ibr-1)*kbl + 1;

          for (jbc=ibr+1; jbc<=nbl; ++jbc) {

             jgl = (jbc-1)*kbl + 1;
//
//     doing the block at ( ibr, jbc )
//
             ijblsk = 0;
             for (p=igl; p<=min(igl+kbl-1,n); ++p) {

                aapp = sva(p);
                if (aapp>Zero) {

                   pSkipped = 0;

                   for (q=jgl; q<=min(jgl+kbl-1,n); ++q) {

                      aaqq = sva(q);
                      if (aaqq>Zero) {
                         aapp0 = aapp;

//
//  .. M x 2 Jacobi SVD ..
//
//     Safe Gram matrix computation
//
                         if (aaqq>=One) {
                            if (aapp>=aaqq) {
                               rotOk = (small*aapp)<=aaqq;
                            } else {
                               rotOk = (small*aaqq)<=aapp;
                            }
                            if (aapp<(big/aaqq)) {
                               aapq = (A(_,p)*A(_,q)*work(p)*work(q)/aaqq)
                                    / aapp;
                            } else {
                               auto _work = work(_(n+1,n+m));

                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, work(p), _work);
                               aapq = _work*A(_,q)*work(q) / aaqq;
                            }
                         } else {
                            if (aapp>=aaqq) {
                               rotOk = aapp<=(aaqq/small);
                            } else {
                               rotOk = aaqq<=(aapp/small);
                            }
                            if (aapp>(small/aaqq)) {
                               aapq = (A(_,p)*A(_,q)*work(p)*work(q)/aaqq)
                                    / aapp;
                            } else {
                               auto _work = work(_(n+1,n+m));

                               _work = A(_,q);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, work(q), _work);
                               aapq = _work*A(_,p)*work(p)/aapp;
                            }
                         }

                         max_aapq = max(max_aapq, abs(aapq));
//
//     TO rotate or NOT to rotate, THAT is the question ...
//
                         if (abs(aapq)>tol) {
                            notRot = 0;
//[RTD]    ROTATED  = ROTATED + 1
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
                                  t = Half/theta;
                                  fastr(3) =  t*work(p) / work(q);
                                  fastr(4) = -t*work(q) / work(p);
                                  blas::rotm(A(_,p), A(_,q), fastr);
                                  if (rhsVec) {
                                     blas::rotm(V(_,p), V(_,q), fastr);
                                  }
                                  sva(q) = aaqq
                                         * sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero,One-t*aqoap*aapq));
                                  max_sinj = max(max_sinj, abs(t));
                               } else {
//
//              .. choose correct signum for THETA and rotate
//
                                  thetaSign = -sign(One,aapq);
                                  if (aaqq>aapp0) {
                                      thetaSign = -thetaSign;
                                  }
                                  t = One
                                    / (theta+thetaSign*sqrt(One+theta*theta));
                                  cs = sqrt(One / (One + t*t));
                                  sn = t*cs;
                                  max_sinj = max(max_sinj, abs(sn));
                                  sva(q) = aaqq
                                         * sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));

                                  apoaq = work(p) / work(q);
                                  aqoap = work(q) / work(p);
                                  if (work(p)>=One) {

                                     if (work(q)>=One) {
                                        fastr(3) =  t*apoaq;
                                        fastr(4) = -t*aqoap;
                                        work(p) *= cs;
                                        work(q) *= cs;
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
                                        work(p) *= cs;
                                        work(q) /= cs;
                                     }
                                  } else {
                                     if (work(q)>=One) {
                                        A(_,q) += t*apoaq*A(_,p);
                                        A(_,p) -= cs*sn*aqoap*A(_,q);
                                        if (rhsVec) {
                                           V(_,q) += t*apoaq*V(_,p);
                                           V(_,p) -= cs*sn*aqoap*V(_,q);
                                        }
                                        work(p) /= cs;
                                        work(q) *= cs;
                                     } else {
                                        if (work(p)>=work(q)) {
                                           A(_,p) -= t*aqoap*A(_,q);
                                           A(_,q) += cs*sn*apoaq*A(_,p);
                                           work(p) *= cs;
                                           work(q) /= cs;
                                           if (rhsVec) {
                                              V(_,p) -= t*aqoap*V(_,q);
                                              V(_,q) += cs*sn*apoaq*V(_,p);
                                           }
                                        } else {
                                           A(_,q) += t*apoaq*A(_,p);
                                           A(_,p) -= cs*sn*aqoap*A(_,q);
                                           work(p) /= cs;
                                           work(q) *= cs;
                                           if (rhsVec) {
                                              V(_,q) += t*apoaq*V(_,p);
                                              V(_,p) -= cs*sn*aqoap*V(_,q);
                                           }
                                        }
                                     }
                                  }
                               }

                            } else {
                               auto _work = work(_(n+1,lWork));

                               if (aapp>aaqq) {
                                  _work = A(_,p);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aapp, One, _work);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aaqq, One, A(_,q));
                                  tmp = -aapq*work(p) / work(q);
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
                                  tmp = -aapq*work(q) / work(p);
                                  A(_,p) += tmp * _work;
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        One, aapp, A(_,p));
                                  sva(p) = aapp*sqrt(max(Zero, One-aapq*aapq));
                                  max_sinj = max(max_sinj, safeMin);
                               }
                            }
//        END IF ROTOK THEN ... ELSE
//
//        In the case of cancellation in updating SVA(q)
//        .. recompute SVA(q)
                            if (pow(sva(q)/aaqq,2)<=rootEps) {
                               if ((aaqq<rootBig) && (aaqq>rootSafeMin)) {
                                  sva(q) = blas::nrm2(A(_,q))*work(q);
                               } else {
                                  t = Zero;
                                  aaqq = One;
                                  lassq(A(_,q), t, aaqq);
                                  sva(q) = t*sqrt(aaqq)*work(q);
                               }
                            }
                            if (pow(aapp/aapp0,2)<=rootEps) {
                               if ((aapp<rootBig) && (aapp>rootSafeMin)) {
                                  aapp = blas::nrm2(A(_,p))*work(p);
                               } else {
                                  t = Zero;
                                  aapp = One;
                                  lassq(A(_,p), t, aapp);
                                  aapp = t*sqrt(aapp)*work(p);
                               }
                               sva(p) = aapp;
                            }
//             end of OK rotation
                         } else {
                            ++notRot;
//[RTD]      SKIPPED  = SKIPPED  + 1
                            ++pSkipped;
                            ++ijblsk;
                         }
                      } else {
                         ++notRot;
                         ++pSkipped;
                         ++ijblsk;
                      }

                      if ((i<=swBand) && (ijblsk>=blSkip)) {
                         sva(p) = aapp;
                         notRot = 0;
                         goto jbcLoopExit;
                      }
                      if ((i<=swBand) && (pSkipped>rowSkip)) {
                         aapp = -aapp;
                         notRot = 0;
                         break;
                      }

                   }
//       end of the q-loop

                   sva(p) = aapp;

                } else {

                   if (aapp==Zero) {
                      notRot += min(jgl+kbl-1,n) -jgl + 1;
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
//  bailed out of the jbc-loop
          for (p=igl; p<=min(igl+kbl-1,n); ++p) {
             sva(p) = abs(sva(p));
          }
//**
       }
//  end of the ibr-loop
//
//  .. update SVA(N)
       if ((sva(n)<rootBig) && (sva(n)>rootSafeMin)) {
          sva(n) = blas::nrm2(A(_,n))*work(n);
       } else {
          t = Zero;
          aapp = One;
          lassq(A(_,n), t, aapp);
          sva(n) = t*sqrt(aapp)*work(n);
       }
//
//  Additional steering devices
//
       if ((i<swBand) && ((max_aapq<=rootTol) || (iswRot<=n))) {
          swBand = i;
       }

       if (i>swBand+1 && max_aapq<sqrt(ElementType(n))*tol
        && ElementType(n)*max_aapq*max_sinj<tol) {
           converged = true;
           break;
       }

       if (notRot>=emptsw) {
           converged = true;
           break;
       }

    }
//  end i=1:NSWEEP loop
//
    if (converged) {
//#:) INFO = 0 confirms successful iterations.
        info = 0;
    } else {
//#:( Reaching this point means that the procedure has not converged.
        info = nSweep - 1;
    }
//
//  Sort the singular values and find how many are above
//  the underflow threshold.
//
    IndexType n2 = 0;
    IndexType n4 = 0;
    for (IndexType p=1; p<=n-1; ++p) {
        const IndexType q = blas::iamax(sva(_(p,n))) + p - 1;
        if (p!=q) {
            swap(sva(p), sva(q));
            swap(work(p), work(q));
            blas::swap(A(_,p), A(_,q));
            if (rhsVec) {
                blas::swap(V(_,p), V(_,q));
            }
        }
        if (sva(p)!=Zero) {
            ++n4;
            if (sva(p)*skl>safeMin) {
                ++n2;
            }
        }
    }
    if (sva(n)!=Zero) {
        ++n4;
        if (sva(n)*skl>safeMin) {
            ++n2;
        }
    }
//
//  Normalize the left singular vectors.
//
    if (lhsVec || controlU) {
        for (IndexType p=1; p<=n2; ++p) {
            A(_,p) *= work(p)/sva(p);
        }
    }
//
//  Scale the product of Jacobi rotations (assemble the fast rotations).
//
    if (rhsVec) {
        if (applyV) {
            for (IndexType p=1; p<=n; ++p) {
                V(_,p) *= work(p);
            }
        } else {
            for (IndexType p=1; p<=n; ++p) {
                V(_,p) *= One / blas::nrm2(V(_,p));
            }
        }
    }
//
//  Undo scaling, if necessary (and possible).
    if (((skl>One) && (sva(1)<big/skl))
     || ((skl<One) && (sva(n2)>safeMin/skl)))
    {
        sva *= skl;
        skl = One;
    }
//
    work(1) = skl;
//  The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
//  then some of the singular values may overflow or underflow and
//  the spectrum is given in this factored representation.
//
    work(2) = ElementType(n4);
//  N4 is the number of computed nonzero singular values of A.
//
    work(3) = ElementType(n2);
//  N2 is the number of singular values of A greater than SFMIN.
//  If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
//  that may carry some information.
//
    work(4) = ElementType(i);
//  i is the index of the last sweep before declaring convergence.
//
    work(5) = max_aapq;
//  MXAAPQ is the largest absolute value of scaled pivots in the
//  last sweep
//
    work(6) = max_sinj;
//  MXSINJ is the largest absolute value of the sines of Jacobi angles
//  in the last sweep
//
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj_impl(SVJ::TypeA                typeA,
         SVJ::JobU                 jobU,
         SVJ::JobV                 jobV,
         GeMatrix<MA>              &A,
         DenseVector<VSVA>         &sva,
         GeMatrix<MV>              &V,
         DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::gesvj<IndexType>(getF77Char(typeA),
                                       getF77Char(jobU),
                                       getF77Char(jobV),
                                       A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       sva.data(),
                                       V.numRows(),
                                       V.data(),
                                       V.leadingDimension(),
                                       work.data(),
                                       work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj_(SVJ::TypeA                typeA,
     SVJ::JobU                 jobU,
     SVJ::JobV                 jobV,
     GeMatrix<MA>              &A,
     DenseVector<VSVA>         &sva,
     GeMatrix<MV>              &V,
     DenseVector<VWORK>        &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(m>=n);

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
        ASSERT(work.length()>=max(IndexType(6),m+n));
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org     = A;
    typename DenseVector<VSVA>::NoView  sva_org   = sva;
    typename GeMatrix<MV>::NoView       V_org     = V;
    typename DenseVector<VWORK>::NoView work_org  = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::svj_impl(typeA, jobU, jobV,
                                             A, sva, V, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic     = A;
    typename DenseVector<VSVA>::NoView  sva_generic   = sva;
    typename GeMatrix<MV>::NoView       V_generic     = V;
    typename DenseVector<VWORK>::NoView work_generic  = work;
//
//  restore output arguments
//
    A    = A_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::svj_impl(typeA, jobU, jobV, A, sva, V, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
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
        std::cerr << "error in: svj.tcc ("
                  << ", m = " << A.numRows()
                  << ", n = " << A.numCols()
                  << ", typeA = " << char(typeA)
                  << ", jobU = " << char(jobU)
                  << ", jobV = " << char(jobV)
                  << ", info = " << info
                  << ") " << std::endl;
        ASSERT(0);
    } else {
        /*
        std::cerr << "passed: svj.tcc ("
                  << ", m = " << A.numRows()
                  << ", n = " << A.numCols()
                  << ", typeA = " << char(typeA)
                  << ", jobU = " << char(jobU)
                  << ", jobV = " << char(jobV)
                  << ", info = " << info
                  << ") " << std::endl;
        */
    }
#   endif

    return info;
}

template <typename MA, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj(SVJ::TypeA                typeA,
    SVJ::JobU                 jobU,
    SVJ::JobV                 jobV,
    GeMatrix<MA>              &A,
    DenseVector<VSVA>         &sva,
    GeMatrix<MV>              &V,
    DenseVector<VWORK>        &work)
{
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org    = A;
    typename DenseVector<VSVA>::NoView  sva_org  = sva;
    typename GeMatrix<MV>::NoView       V_org    = V;
    typename DenseVector<VSVA>::NoView  work_org = work;

    svj_(SVJ::Lower, jobU, jobV, A, sva, V, work);

    A    = A_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;

    svj_(SVJ::Upper, jobU, jobV, A, sva, V, work);

    A    = A_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;
    svj_(SVJ::General, jobU, jobV, A, sva, V, work);

    A    = A_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;
#   endif

    return svj_(typeA, jobU, jobV, A, sva, V, work);
}


template <typename MA, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj(SVJ::JobU                 jobU,
    SVJ::JobV                 jobV,
    GeMatrix<MA>              &A,
    DenseVector<VSVA>         &sva,
    GeMatrix<MV>              &V,
    DenseVector<VWORK>        &work)
{
    return svj(SVJ::General, jobU, jobV, A, sva, V, work);
}

template <typename MA, typename VSVA, typename MV, typename VWORK>
typename TrMatrix<MA>::IndexType
svj(SVJ::JobU                 jobU,
    SVJ::JobV                 jobV,
    TrMatrix<MA>              &A,
    DenseVector<VSVA>         &sva,
    GeMatrix<MV>              &V,
    DenseVector<VWORK>        &work)
{
    SVJ::TypeA  upLo = (A.upLo()==Upper) ? SVJ::Upper : SVJ::Lower;

    return svj(upLo, jobU, jobV, A.general(), sva, V, work);
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VSVA, typename MV, typename VWORK>
typename MA::IndexType
svj(SVJ::TypeA                typeA,
    SVJ::JobU                 jobU,
    SVJ::JobV                 jobV,
    MA                        &&A,
    VSVA                      &&sva,
    MV                        &&V,
    VWORK                     &&work)
{
    typename MA::IndexType info;

    CHECKPOINT_ENTER;
    info = svj(typeA, jobU, jobV, A, sva, V, work);
    CHECKPOINT_LEAVE;

    return info;
}

template <typename MA, typename VSVA, typename MV, typename VWORK>
typename MA::IndexType
svj(SVJ::JobU                 jobU,
    SVJ::JobV                 jobV,
    MA                        &&A,
    VSVA                      &&sva,
    MV                        &&V,
    VWORK                     &&work)
{
    typename MA::IndexType info;

    CHECKPOINT_ENTER;
    info = svj(jobU, jobV, A, sva, V, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_SVJ_TCC
