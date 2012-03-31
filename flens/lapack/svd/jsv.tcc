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
       SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,
      $                   M, N, A, LDA, SVA, U, LDU, V, LDV,
      $                   WORK, LWORK, IWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                    --
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

#ifndef FLENS_LAPACK_SVD_JSV_TCC
#define FLENS_LAPACK_SVD_JSV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
/*
template <typename MA, typename VSVA, typename MU, typename MV,
          typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
jsv_generic(JSV::Accuracy             accuracy,
            JSV::JobU                 jobU,
            JSV::JobV                 jobV,
            bool                      restrictedRange,
            bool                      considerTransA,
            bool                      perturb,
            GeMatrix<MA>              &A,
            DenseVector<VSVA>         &sva,
            GeMatrix<MU>              &U,
            GeMatrix<MV>              &V,
            DenseVector<VWORK>        &work,
            DenseVector<VIWORK>       &iwork)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType Zero(0), One(1);

    const Underscore<IndexType>  _;

    const bool lsvec  = (jobU==ComputeU) || (jobU==FullsetU);
    const bool jracc  = (jobV=='J');
    const bool rsvec  = (jobV=='V') || jracc;
    const bool rowpiv = (jobA=='F') || (jobA=='G');
    const bool l2rank = (jobA=='R');
    const bool l2aber = (jobA=='A');
    const bool errest = (jobA=='E') || (jobA=='G');
    const bool l2tran = (jobT=='T');
    const bool l2kill = (jobR=='R');
    const bool defr   = (jobR=='N');
    const bool l2pert = (jobP=='P');

    IndexType info = 0;

//
//  Quick return for void matrix (Y3K safe)
//#:)
    if (m==0 || n==0) {
        return info;
    }
//
//  Set numerical parameters
//
//! NOTE: Make sure DLAMCH() does not fail on the target architecture.
//
    const ElementType eps     = lamch<ElementType>(Eps);
    const ElementType safeMin = lamch<ElementType>(SafeMin);
    const ElementType small   = safeMin / eps;
    const ElementType big     = lamch<ElementType>(OverflowThreshold);
//
//  Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
//
//! If necessary, scale SVA() to protect the largest norm from
//  overflow. It is possible that this scaling pushes the smallest
//  column norm left from the underflow threshold (extreme case).
//
    ElementType scaleM = One / sqrt(ElementType(m)*ElementType(n));
    bool noScale = true;
    bool goScale = true;
    for (IndexType p=1; p<=n; ++p) {
        aapp = Zero;
        aaqq = One;
        lassq(A(_,p), aapp, aaqq);
        if (aapp>big) {
            return -9;
        }
        aaqq = sqrt(aaqq);
        if (aapp<(big/aaqq) && noscal) {
            sva(p)  = aapp * aaqq;
        } else {
            noscal  = false;
            sva(p)  = aapp * ( aaqq * scalem )
            if (goScale) {
                goScale = false;
                sva(_(1,p)) *= scaleM;
            }
        }
    }

    if (noScale) {
        scaleM = One;
    }

    aapp = Zero;
    aaqq = big;
    for (IndexType p=1; p<=n) {
        aapp = max(aapp, sva(p));
        if (sva(p)!=Zero) {
            aaqq = min(aaqq, sva(p));
        }
    }
//
//    Quick return for zero M x N matrix
//#:)
    if (aapp==Zero) {
        if (lsvec) {
            U = Zero;
            U.diag(0) = One;
        }
        if (rsvec) {
            V = Zero;
            V.diag(0) = One;
        }
        work(1) = One;
        work(2) = One;
        if (errest) {
            work(3) = One;
        }
        if (lsvec && rsvec) {
            work(4) = One;
            work(5) = One;
        }
        if (l2tran) {
            work(6) = Zero;
            work(7) = Zero;
        }
        iwork(1) = 0;
        iwork(2) = 0;
        iwork(3) = 0;
        return info;
    }
//
//  Issue warning if denormalized column norms detected. Override the
//  high relative accuracy request. Issue licence to kill columns
//  (set them to zero) whose norm is less than sigma_max / BIG (roughly).
//#:(
    warning = 0
    if (aaqq<=SFMIN) {
        l2rank = true;
        l2kill = true;
        warning = 1;
    }
//
//  Quick return for one-column matrix
//#:)
    auto U1 = U(_,_(1,n));
    if (n==1) {

        if (lsvec) {
            lascl(LASCL::FullMatrix, 0, 0, sva(1), scaleM, A);
            U1 = A;
//          computing all M left singular vectors of the M x 1 matrix
            if (nu!=n) {
                auto tau   = work(_(1,n));
                auto _work = work(_(n+1,lWork));
                qrf(U1, tau, _work);
                orgqr(1, U, tau, _work);
                U1 = A;
            }
        }
        if (rsvec) {
            V(1,1) = One
        }
        if (sva(1)<big*scaleM) {
            sva(1) /= scaleM;
            scalem = One;
        }
        work(1) = One / scaleM;
        work(2) = One;
        if (sva(1)!=Zero) {
            iwork(1) = 1;
            if (sva(1)/scaleM>=safeMin) {
                iwork(2) = 1;
            } else {
                iwork(2) = 0;
            }
        } else {
            iwork(1) = 0;
            iwork(2) = 0;
        }
        if (errest) {
            work(3) = One;
        }
        if (lsvec && rsvec) {
            work(4) = One;
            work(5) = One;
        }
        if (l2tran) {
            work(6) = Zero;
            work(7) = Zero;
        }
        return info;

    }

    bool transp = false;
    bool l2tran = l2tran && m==n;

    ElementType aatMax = -One;
    ElementType aatMin = big;
    if (rowpiv || l2tran) {
//
//      Compute the row norms, needed to determine row pivoting sequence
//      (in the case of heavily row weighted A, row pivoting is strongly
//      advised) and to collect information needed to compare the
//      structures of A * A^t and A^t * A (in the case L2TRAN.EQ.true).
//
        if (l2tran) {
            for (IndexType p=1; p<=m; ++p) {
                xsc = Zero;
                tmp = One;
                lassq(A(p,_), xsc, tmp);
//              DLASSQ gets both the ell_2 and the ell_infinity norm
//              in one pass through the vector
                work(m+n+p) = xsc * scaleM;
                work(n+p)   = xsc * (scaleN*sqrt(tmp));
                aatMax      = max(aatMax, work(n+p));
                if (work(n+p)!=Zero) {
                    aatMin = min(aatMin, work(n+p));
                }
            }
        } else {
            for (IndexType p=1; p<=m; ++p) {
                const IndexType jp = blas::iamax(A(p,_));
                work(m+n+p) = scaleM*abs(A(p,jp));
                aatMax = max(aatMax, work(m+n+p));
                aatMin = min(aatMin, work(m+n+p));
            }
        }

    }
//
//  For square matrix A try to determine whether A^t  would be  better
//  input for the preconditioned Jacobi SVD, with faster convergence.
//  The decision is based on an O(N) function of the vector of column
//  and row norms of A, based on the Shannon entropy. This should give
//  the right choice in most cases when the difference actually matters.
//  It may fail and pick the slower converging side.
//
    entra  = Zero
    entrat = Zero
    if (l2tran) {

        xsc = Zero;
        tmp = One;
        lassq(sva, xsc, tmp);
        tmp = One / tmp;

        entra = Zero
        for (IndexType p=1; p<=n; ++p) {
            const ElementType _big = pow(sva(p)/xsc, 2) * tmp;
            if (_big!=Zero) {
                entra += _big * log(_big);
            }
        }
        entra = - entra / log(ElementType(n));
//
//      Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
//      It is derived from the diagonal of  A^t * A.  Do the same with the
//      diagonal of A * A^t, compute the entropy of the corresponding
//      probability distribution. Note that A * A^t and A^t * A have the
//      same trace.
//
        entrat = Zero
        for (IndexType p=n+1; p<=n+m; ++p) {
            const ElementType _big = pow(work(p)/xsc, 2) * tmp;
            if (_big!=Zero) {
                entrat += _big * log(_big);
            }
        }
        entrat = -entrat / log(ElementType(m));
//
//      Analyze the entropies and decide A or A^t. Smaller entropy
//      usually means better input for the algorithm.
//
        transp = entrat<entra;
//
//      If A^t is better than A, transpose A.
//
        if (transp) {
//          In an optimal implementation, this trivial transpose
//          should be replaced with faster transpose.
//          TODO: in-place transpose:
//                transpose(A);
            for (IndexType p=1; p<=n-1; ++p) {
                for (IndexType q=p+1; q<=n; ++q) {
                    swap(A(q,p), A(p,q));
                }
            }
            for (IndexType p=1; p<=n; ++p) {
                work(m+n+p) = sva(p);
                sva(p)      = work(n+p);
            }
            swap(aapp, aatMax);
            swap(aaqq, aatMin);
            swap(lsvec, rsvec);
            if (lsvec) {
                // Lehn: transposing A is only considered if A is square,
                //       so in this case m==n and therefore nu==n.  Or am I
                //       wrong??
                ASSERT(nu==n);
            }

            rowpiv = true;
        }

    }
//
//  Scale the matrix so that its maximal singular value remains less
//  than DSQRT(BIG) -- the matrix is scaled so that its maximal column
//  has Euclidean norm equal to DSQRT(BIG/N). The only reason to keep
//  DSQRT(BIG) instead of BIG is the fact that DGEJSV uses LAPACK and
//  BLAS routines that, in some implementations, are not capable of
//  working in the full interval [SFMIN,BIG] and that they may provoke
//  overflows in the intermediate results. If the singular values spread
//  from SFMIN to BIG, then DGESVJ will compute them. So, in that case,
//  one should use DGESVJ instead of DGEJSV.
//
    const ElementType bigRoot = sqrt(big);
    tmp = sqrt(big/ElementType(n));

    lascl(LASCL::FullMatrix, 0, 0, aapp, tmp, sva);
    if (aaqq>aapp*safeMin) {
        aaqq = (aaqq/aapp) * tmp;
    } else {
        aaqq = (aaqq*tmp) / aapp;
    }
    tmp *= scaleM;
    lascl(LASCL::FullMatrix, 0, 0, aapp, tmp, A);
//
//  To undo scaling at the end of this procedure, multiply the
//  computed singular values with USCAL2 / USCAL1.
//
    uScale1 = tmp;
    uScale2 = aapp;

    if (l2kill) {
//      L2KILL enforces computation of nonzero singular values in
//      the restricted range of condition number of the initial A,
//      sigma_max(A) / sigma_min(A) approx. DSQRT(BIG)/DSQRT(SFMIN).
        xsc = sqrt(safeMin);
    } else {
        xsc = small;
//
//      Now, if the condition number of A is too big,
//      sigma_max(A) / sigma_min(A)>DSQRT(BIG/N) * EPSLN / SFMIN,
//      as a precaution measure, the full SVD is computed using DGESVJ
//      with accumulated Jacobi rotations. This provides numerically
//      more robust computation, at the cost of slightly increased run
//      time. Depending on the concrete implementation of BLAS and LAPACK
//      (i.e. how they behave in presence of extreme ill-conditioning) the
//      implementor may decide to remove this switch.
        if (aaqq<sqrt(safeMin) && lsvec && rsvec ) {
            jracc = true;
        }

    }
    if (aaqq<xsc) {
        for (IndexType p=1; p<=n; ++p) {
            if (sva(p)<xsc) {
                A(_,p) = Zero;
                sva(p) = Zero;
            }
        }
    }
//
//  Preconditioning using QR factorization with pivoting
//
    if (rowpiv) {
//      Optional row permutation (Bjoerck row pivoting):
//      A result by Cox and Higham shows that the Bjoerck's
//      row pivoting combined with standard column pivoting
//      has similar effect as Powell-Reid complete pivoting.
//      The ell-infinity norms of A are made nonincreasing.
        for (IndexType p=1; p<=m-1; ++p) {
            const IndexType q = blas::iamax(work(_(m+n+p,2*m+n))) + p - 1;
            iwork(2*n+p) = q;
            if (p!=q) {
                swap(work(m+n+p), work(m+n+q));
            }
        }
        laswp(A, iwork(_(2*n+1, 2*n+m-1)));
    }
//
//  End of the preparation phase (scaling, optional sorting and
//  transposing, optional flushing of small columns).
//
//  Preconditioning
//
//  If the full SVD is needed, the right singular vectors are computed
//  from a matrix equation, and for that we need theoretical analysis
//  of the Businger-Golub pivoting. So we use DGEQP3 as the first RR QRF.
//  In all other cases the first RR QRF can be chosen by other criteria
//  (eg speed by replacing global with restricted window pivoting, such
//  as in SGEQPX from TOMS # 782). Good results will be obtained using
//  SGEQPX with properly (!) chosen numerical parameters.
//  Any improvement of DGEQP3 improves overal performance of DGEJSV.
//
//  A * P1 = Q1 * [ R1^t 0]^t:
    auto _tau   = work(_(1,n));
    auto _work  = work(_(n+1,lWork));
    auto _iwork = iwork(_(1,n));

//      .. all columns are free columns
    _iwork = 0;

    qp3(A, _iwork, _tau, _work);
//
//  The upper triangular matrix R1 from the first QRF is inspected for
//  rank deficiency and possibilities for deflation, or possible
//  ill-conditioning. Depending on the user specified flag L2RANK,
//  the procedure explores possibilities to reduce the numerical
//  rank by inspecting the computed upper triangular factor. If
//  L2RANK or l2aber are up, then DGEJSV will compute the SVD of
//  A + dA, where ||dA|| <= f(M,N)*EPSLN.
//
    IndexType nr = 1;
    if (l2aber) {
//      Standard absolute error bound suffices. All sigma_i with
//      sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
//      agressive enforcement of lower numerical rank by introducing a
//      backward error of the order of N*EPSLN*||A||.
        tmp = sqrt(ElementType(n))*eps;
        for (IndexType p=2; p<=n; ++p) {
            if (abs(A(p,p))>=tmp*abs(A(1,1))) {
                ++nr;
            } else {
                break;
            }
        }
    } else if (l2rank) {
//      .. similarly as above, only slightly more gentle (less agressive).
//      Sudden drop on the diagonal of R1 is used as the criterion for
//      close-to-rank-defficient.
        tmp = sqrt(safeMin);
        for (IndexType p=2; p<=n; ++p) {
            if (abs(A(p,p))<eps*abs(A(p-1,p-1)) || abs(A(p,p))<small
             || (l2kill && abs(A(p,p))<tmp))
            {
                break;
            }
            ++nr;
        }

    } else {
//      The goal is high relative accuracy. However, if the matrix
//      has high scaled condition number the relative accuracy is in
//      general not feasible. Later on, a condition number estimator
//      will be deployed to estimate the scaled condition number.
//      Here we just remove the underflowed part of the triangular
//      factor. This prevents the situation in which the code is
//      working hard to get the accuracy not warranted by the data.
        tmp = sqrt(safeMin);
        for (IndexType p=2; p<=n; ++p) {
            if (abs(A(p,p))<small || (l2kill && abs(A(p,p))<tmp)) {
                break;
            }
            ++nr;
        }

    }

    bool almort = false;
    if (nr==n) {
        ElementType maxprj = One;
        for (IndexType p=2; p<=n; ++p) {
            maxprj = min(maxprj, abs(A(p,p))/sva(iwork(p)));
        }
        if (pow(maxprj,2)>=One-ElementType(n)*eps) {
            almort = true;
        }
    }


    sconda = -One;
    condr1 = -One;
    condr2 = -One;

    if (errest) {
        if (n==nr) {
            if (rsvec) {
//              .. V is available as workspace
                V.upper() = A(_(1,n),_).upper();
                for (IndexType p=1; p<=n; ++p) {
                    tmp = sva(iwork(p));
                    V(_(1,p),p) *= One/tmp;
                }
                auto _work  = work(_(n+1, 4*n));
                auto _iwork = iwork(_(2*n+m+1, 3*n+m));
                pocon(V.upper(), One, tmp, _work, _iwork);
            } else if (lsvec) {
//              .. U is available as workspace
                auto _U = U(_(1,n),_(1,n));
                _U.upper() = A(_(1,n),_).upper();
                for (IndexType p=1; p<=n; ++p) {
                    tmp = sva(iwork(p));
                    U(_(1,p),p) *= One/tmp;
                }
                auto _work  = work(_(n+1, 4*n));
                auto _iwork = iwork(_(2*n+m+1, 3*n+m));
                pocon(_U.upper(), One, tmp, _work, _iwork);
            } else {
                auto                       _work1 = work(_(n+1,n+n*n));
                auto                       _work2 = work(_(n+n*n+1,n+n*n+3*n));
                auto                       _iwork = work(_(2*n+m+1,3*n+m));
                GeMatrixView<ElementType>  Work(n, n, _work1, n);

                Work.upper() = A(_(1,n),_).upper();
                for (IndexType p=1; p<=n; ++p) {
                    tmp = sva(iwork(p));
                    Work(_(1,p),p) *= One/tmp;
                }
//              .. the columns of R are scaled to have unit Euclidean lengths.
                pocon(Work.upper(), One, tmp, _work2, _iwork);
            }
            sconda = One/sqrt(tmp);
//          SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).
//          N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
        } else {
            sconda = -One;
        }
    }

    l2pert = l2pert && abs(A(1,1)/A(nr,nr))>sqrt(bigRoot);
//  If there is no violent scaling, artificial perturbation is not needed.
//
//  Phase 3:
//
    if (!(rsvec || lsvec)) {
//
//      Singular Values only
//
//      .. transpose A(1:NR,1:N)
        for (IndexType p=1; p<=min(n-1, nr); ++p) {
            A(p,_(p+1,n)) = A(_(p+1,n),p);
        }
//
//      The following two DO-loops introduce small relative perturbation
//      into the strict upper triangle of the lower triangular matrix.
//      Small entries below the main diagonal are also changed.
//      This modification is useful if the computing environment does not
//      provide/allow FLUSH TO Zero underflow, for it prevents many
//      annoying denormalized numbers in case of strongly scaled matrices.
//      The perturbation is structured so that it does not introduce any
//      new perturbation of the singular values, and it does not destroy
//      the job done by the preconditioner.
//      The licence for this perturbation is in the variable L2PERT, which
//      should be false if FLUSH TO Zero underflow is active.
//
        if (! almort) {

            if (l2pert) {
//              XSC = DSQRT(SMALL)
                xsc = eps / ElementType(n);
                for (IndexType q=1; q<=nr; ++q) {
                    tmp = xsc*abs(A(q,q));
                    for (IndexType p=1; p<=n, ++p) {
                        if ((p>q && abs(A(p,q))<=tmp) || p<q) {
                            A(p,q) = sign(tmp, A(p,q));
                        }
                    }
                }
            } else {
                A(_(1,nr),_(1,nr)).strictUpper() = Zero;
            }
//
//          .. second preconditioning using the QR factorization
//
            auto _A    = A(_(1,n),_(1,nr));
            auto _tau  = work(_(1,nr));
            auto _work = work(_(n+1,lWork));
            qrf(_A, _tau, _work);
//
//          .. and transpose upper to lower triangular
            for (IndexType p=1; p<=nr-1; ++p) {
                A(p,_(p+1,nr)) = A(_(p+1,nr),p);
            }

        }
//
//      Row-cyclic Jacobi SVD algorithm with column pivoting
//
//      .. again some perturbation (a "background noise") is added
//      to drown denormals
        if (l2pert) {
//          XSC = DSQRT(SMALL)
            xsc = eps / ElementType(n);
            for (IndexType q=1; q<=nr; ++q) {
                ElementType tmp = xsc*abs(A(q,q));
                for (IndexType p=1; p<=nr; ++p) {
                    if (((p>q) && abs(A(p,q))<=tmp) || p<q) {
                        A(p,q) = sign(tmp, A(p,q));
                    }
                }
            }
        } else {
            A(_(1,nr),_(1,nr)).strictUpper() = Zero;
        }
//
//      .. and one-sided Jacobi rotations are started on a lower
//      triangular matrix (plus perturbation which is ignored in
//      the part which destroys triangular form (confusing?!))
//
        auto _A    = A(_(1,nr),_(1,nr));
        auto _sva  = sva(_(1,nr));

        svj(SVJ::Lower, SVJ::NoU, SVJ::NoV, _A, _sva, V, work);

        scaleM  = work(1);
        numRank = nint(work(2));


    } else if (rsvec && !lsvec) {
//
//      -> Singular Values and Right Singular Vectors <-
//
        if (almort) {
//
//          .. in this case NR equals N
            ASSERT(nr==n);
            for (IndexType p=1; p<=nr; ++p) {
                V(_(p,n),p) = A(p,_(p,n));
            }

            auto _V = V(_,_(1,nr));
            _V.strictUpper() = Zero;
            svj(SVJ::Lower, SVJ::ComputeU, SVJ::NoV, _V, sva, A, work);

            scaleM  = work(1);
            numRank = nint(work(2));

        } else {
//
//      .. two more QR factorizations ( one QRF is not enough, two require
//      accumulated product of Jacobi rotations, three are perfect )
//
            A(_(1,nr),_(1,nr)).strictLower() = Zero;

            auto _A     = A(_(1,nr),_);
            auto _tau1  = work(_(1,nr));
            auto _work1 = work(_(n+1,lWork));
            lqf(_A, _tau1, _work1);

            auto _V = V(_(1,nr),_(1,nr))
            _V.lower()       = _A(_,_(1,nr)).lower();
            _V.strictUpper() = Zero;

            auto _tau2  = work(_(n+1,n+nr));
            auto _work2 = work(_(2*n+1,lWork));
            qrf(_V, _tau2, _work2);

            for (IndexType p=1; p<=nr; ++p) {
                V(_(p,nr),p) = V(p,_(p,nr));
            }
            _V.strictUpper() = Zero;

            auto _sva  = sva(_(1,nr));
            svj(SVJ::Lower, SVJ::ComputeU, SVJ::NoV, _V, _sva, U, _work1);

            scaleM  = work(n+1);
            numRank = nint(work(n+2));
            if (nr<n) {
                V(_(nr+1,n),_(1,nr)) = Zero;
                V(_(1,nr),_(nr+1,n)) = Zero;

                V(_(nr+1,n),_(nr+1,n)) = Zero;
                V(_(nr+1,n),_(nr+1,n)).diag(0) = One;
            }

            ormlq(Left, Trans, _A, _tau1, V, _work1);

        }

        for (IndexType p=1; p<=n; ++p) {
            A(iwork(p),_) = V(p,_);
        }
        V = A;

        if (transp) {
            U = V;
        }

    } else if (lsvec && !rsvec) {
//
//      .. Singular Values and Left Singular Vectors                 ..
//
//      .. second preconditioning step to avoid need to accumulate
//      Jacobi rotations in the Jacobi iterations.

        auto _U = U(_(1,nr),_(1,nr));
        auto _tau   = work(_(n+1,n+nr));
        auto _sva   = sva(_(1,nr));
        auto _work1 = work(_(n+1,lWork));
        auto _work2 = work(_(2*n+1,lWork));

        for (IndexType p=1; p<=nr; ++p) {
            A(p,_(p,nr)) = U(_(p,nr),p);
        }
        _U.strictUpper() = Zero;

        qrf(U(_(1,n),_(1,nr)), _tau, _work2);

        for (IndexType p=1; p<=nr-1; ++p) {
            U(p,_(p+1,nr)) = U(_(p+1,nr),p);
        }
        _U.strictUpper() = Zero;

        svj(SVJ::Lower, SVJ::ComputeU, SVJ::NoV, _U, _sva, A, _work1);
        scaleM  = work(n+1);
        numRank = nint(work(n+2));

        if (nr<m) {
            U(_(nr+1,m),_(1,nr)) = Zero;
            if (nr<nu) {
                U(_(1,nr),_(nr+1,nu))           = Zero;

                U(_(nr+1,m),_(nr+1,nu))         = Zero;
                U(_(nr+1,m),_(nr+1,nu)).diag(0) = One;
            }
        }

        ormqr(Left, NoTrans, A, work, U);

        if (rowpiv) {
            auto piv = iwork(_(2*n+1,2*n+m-1));
            laswp(U, piv.reverse());
        }

        for (IndexType p=1; p<=nu; ++p) {
            xsc = One / blas::nrm2(U(_,p));
            U(_,p) *= xsc;
        }

        if (transp) {
            V = U;
        }
//
    } else {
//
//      .. Full SVD ..
//
        if (!jracc) {

            if (!almort) {
//
//              Second Preconditioning Step (QRF [with pivoting])
//              Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
//              equivalent to an LQF CALL. Since in many libraries the QRF
//              seems to be better optimized than the LQF, we do explicit
//              transpose and use the QRF. This is subject to changes in an
//              optimized implementation of DGEJSV.
//
                for (IndexType p=1; p<=nr; ++p) {
                    V(_(p,n),p) = A(p,_(p,n));
                }
//
//              .. the following two loops perturb small entries to avoid
//              denormals in the second QR factorization, where they are
//              as good as zeros. This is done to avoid painfully slow
//              computation with denormals. The relative size of the
//              perturbation is a parameter that can be changed by the
//              implementer. This perturbation device will be obsolete on
//              machines with properly implemented arithmetic.
//              To switch it off, set L2PERT=false To remove it from  the
//              code, remove the action under L2PERT=true, leave the ELSE
//              part. The following two loops should be blocked and fused with
//              the transposed copy above.
//
                if (l2pert) {
                    xsc = sqrt(small);
                    for (IndexType q=1; q<=nr; ++q) {
                        tmp = xsc*abs(V(q,q));
                        for (IndexType p=1; p<=n; ++p) {
                            if (p>q && abs(V(p,q))<=tmp || p<q) {
                                V(p,q) = sign(tmp, V(p,q));
                            }
                            if (p<q) {
                                V(p,q) = -V(p,q);
                            }
                        }
                    }
                } else {
                    V(_(1,nr),_(1,nr)).strictUpper() = Zero;
                }
//
//              Estimate the row scaled condition number of R1
//              (If R1 is rectangular, N > NR, then the condition number
//              of the leading NR x NR submatrix is estimated.)
//
                auto          _work1 = work(_(2*n+1,2*n+nr*nr));
                auto          _work2 = work(_(2*n+nr*nr+1, 2*n+nr*nr+3*nr));
                auto          _iwork = iwork(_(m+2*n+1,m+2*n+nr));
                GeMatrixView<ElementType>  Work(nr, nr, _work1, nr);

                Work.lower() = V.lower();
                for (IndexType p=1; p<=nr; ++p) {
                    tmp = blas::nrm2(Work(_(p,nr),p));
                    Work(_(p,nr),p) *= One/tmp;
                }
                pocon(Work, One, tmp, _work2, _iwork);
                condr1 = One / sqrt(tmp);
//              .. here need a second oppinion on the condition number
//              .. then assume worst case scenario
//              R1 is OK for inverse <=> condr1 < DBLE(N)
//              more conservative    <=> condr1 < DSQRT(DBLE(N))
//
                cond_ok = sqrt(ElementType(nr));
//[TP]          COND_OK is a tuning parameter.

                if (condr1<cond_ok) {
//                  .. the second QRF without pivoting. Note: in an optimized
//                  implementation, this QRF should be implemented as the QRF
//                  of a lower triangular matrix.
//                  R1^t = Q2 * R2
                    auto tau    = work(_(n+1,n+nr));
                    auto _work  = work(_(2*n+1,lWork));

                    qrf(V(_,_(1,nr)), tau, _work);

                    if (l2pert) {
                        xsc = sqrt(small) / Eps;
                        for (IndexType p=2; p<=nr; ++p) {
                            for (IndexType q=1; q<=p-1; ++q) {
                                tmp = xsc*min(abs(V(p,p)), abs(V(q,q)));
                                if (abs(V(q,p))<=tmp) {
                                    V(q,p) = sign(tmp, V(q,p));
                                }
                            }
                        }
                    }
//
                    if (nr!=n) {
                        auto _work = work(_(2*n+1, 2*n+n*nr));
                        GeMatrixView<ElementType>  Work(n, nr, _work, n);

                        Work = V(_,_(1,nr));
                    }
//                  .. save ...
//
//               .. this transposed copy should be better than naive
//                  TODO:  auto _V = V(_(1,nr),_(1,nr));
//                         _V.lower() = transpose(_V.upper());
//
                    for (IndexType p=1; p<=nr-1; ++p) {
                        V(_(p+1,nr),p) = V(p,_(p+1,nr));
                    }

                    condr2 = condr1;

                } else {
//
//                  .. ill-conditioned case: second QRF with pivoting
//                  Note that windowed pivoting would be equaly good
//                  numerically, and more run-time efficient. So, in
//                  an optimal implementation, the next call to DGEQP3
//                  should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
//                  with properly (carefully) chosen parameters.
//
//                  R1^t * P2 = Q2 * R2
                    auto _V    = V(_,_(1,nr));
                    auto piv   = iwork(_(n+1,n+nr));
                    auto tau   = work(_(n+1,,n+nr));
                    auto _work = work(_(2*n+1, lWork));

                    piv = 0;
                    qp3(_V, piv, tau, _work);
                    if (l2pert) {
                        xsc = sqrt(small);
                        for (IndexType p=2; p<=nr; ++p) {
                            for (IndexType q=1; q<=p-1; ++q) {
                                tmp = xsc*min(abs(V(p,p)), abs(V(q,q)));
                                if (abs(V(q,p))<=tmp) {
                                    V(q,p) = sign(tmp, V(q,p));
                                }
                            }
                        }
                    }

                    auto  _work1 = work(_(2*n+1, 2*n+n*nr));
                    GeMatrixView<ElementType>  Work1(n, nr, _work1, n);

                    Work1 = V(_,_(1,nr));

                    if (l2pert) {
                        xsc = sqrt(small);
                        for (IndexType p=2; p<=nr; ++p) {
                            for (IndexType q=1; q<=p-1; ++q) {
                                tmp = xsc * min(abs(V(p,p)), abs(V(q,q)));
                                V(p,q) = -sign(tmp, V(q,p));
                            }
                        }
                    } else {
                        V(_(1,nr),_(1,nr)).strictLower() = Zero;
                    }
//                  Now, compute R2 = L3 * Q3, the LQ factorization.
                    auto _V     = V(_(1,nr),_(1,nr));
                    auto tau    = work(_(2*n+n*nr+1,2*n+n*nr+nr));
                    auto _work2 = work(_(2*n+n*nr+nr+1,lWork));

                    lqf(_V, tau, _work2);
//                  .. and estimate the condition number
                    auto  _work3 = work(_(2*n+n*nr+nr+1,2*n+n*nr+nr+nr*nr));
                    GeMatrixView<ElementType>  Work3(nr, nr, _work3, nr);

                    Work3.lower() = _V.lower();

                    for (IndexType p=1; p<=nr; ++p) {
                        tmp = blas::nrm2(Work3(p,_(1,p)));
                        Work3(p,_(1,p)) *= One/tmp;
                    }
                    auto _work4 = work(_(2*n+n*nr+nr+nr*nr+1,
                                         2*n+n*nr+nr+nr*nr+3*nr));
                    auto _iwork = iwork(_(m+2*n+1, m+2*n+nr));
                    pocon(Work3.lower(), One, tmp, _work4, _iwork);
                    condr2 = One / sqrt(tmp);

                    if (condr2>=cond_ok) {
//                      .. save the Householder vectors used for Q3
//                      (this overwrittes the copy of R2, as it will not be
//                      needed in this branch, but it does not overwritte the
//                      Huseholder vectors of Q2.).
                        Work1(_(1,nr),_(1,nr)).upper() = V.upper();
//                      .. and the rest of the information on Q3 is in
//                      WORK(2*N+N*NR+1:2*N+N*NR+N)
                    }

                }

                if (l2pert) {
                    xsc = sqrt(small);
                    for (IndexType q=2; q<=nr; ++q) {
                        tmp = xsc * V(q,q);
                        for (IndexType p=1; p<=q-1; ++p) {
//                          V(p,q) = - DSIGN( TEMP1, V(q,p) )
                            V(p,q) = -sign(tmp, V(p,q));
                        }
                    }
                } else {
                    V(_(1,nr),_(1,nr)).strictLower();
                }
//
//              Second preconditioning finished; continue with Jacobi SVD
//              The input matrix is lower triangular.
//
//              Recover the right singular vectors as solution of a well
//              conditioned triangular matrix equation.
//
                if (condr1<cond_ok) {
                    auto _U    = U(_(1,nr),_(1,nr));
                    auto _V    = V(_(1,nr),_(1,nr));
                    auto _sva  = sva(_(1,nr));
                    auto _work = work(_(2*n+n*nr+nr+1,lWork));

                    svj(SVJ::Lower, SVJ::ComputeU, SVJ::NoV,
                        _V, _sva, U, _work);
                    scaleM = _work(1);
                    numRank = nint(_work(2));

                    for (IndexType p=1; p<=nr; ++p) {
                        _U(_,p) = _V(_,p);
                        _V(_,p) *= sva(p);
                    }

//                  .. pick the right matrix equation and solve it
//
                    if (nr==n) {
//:))                   .. best case, R1 is inverted. The solution of this
//                      matrix equation is Q2*V2 = the product of the Jacobi
//                      rotations used in DGESVJ, premultiplied with the
//                      orthogonal matrix  from the second QR factorization.
                        const auto _A = A(_(1,nr),_(1,nr));
                        blas::sm(Left, NoTrans, One, _A.upper(), _V);
                    } else {
//                      .. R1 is well conditioned, but non-square. Transpose(R2)
//                      is inverted to get the product of the Jacobi rotations
//                      used in DGESVJ. The Q-factor from the second QR
//                      factorization is then built in explicitly.

                        auto  _work = work(_(2*n+1, 2*n+n*nr));
                        GeMatrixView<ElementType>  Work1(nr, nr, _work, n);
                        GeMatrixView<ElementType>  Work(n, nr, _work, n);

                        blas::sm(Left, NoTrans, One, Work1.upper(), _V);
                        if (nr<n) {
                            V(_(nr+1,n),_(1,nr)) = Zero;
                            V(_(1,nr),_(nr+1,n)) = Zero;
                            V(_(nr+1,n),_(nr+1,n)) = Zero;
                            V(_(nr+1,n),_(nr+1,n)).diag(0) = One;
                        }
                        auto tau = work(_(n+1,n+nr));
                        auto _work_ormqr = work(_(2*n+n*nr+nr+1,lWork));
                        ormqr(Left, NoTrans, Work, tau, V, _work_ormqr);
                    }
//
                } else if (condr2<cond_ok) {
//
//:)                .. the input matrix A is very likely a relative of
//                  the Kahan matrix :)
//                  The matrix R2 is inverted. The solution of the matrix
//                  equation is Q3^T*V3 = the product of the Jacobi rotations
//                  (appplied to the lower triangular L3 from the LQ
//                  factorization of R2=L3*Q3), pre-multiplied with the
//                  transposed Q3.
                    CALL DGESVJ( 'L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U,
     $                  LDU, work(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO )
                    SCALEM  = work(2*N+N*NR+NR+1)
                    numRank = IDNINT(work(2*N+N*NR+NR+2))
                    DO 3870 p = 1, NR
                        CALL DCOPY( NR, V(1,p), 1, U(1,p), 1 )
                        CALL DSCAL( NR, SVA(p),    U(1,p), 1 )
 3870               CONTINUE
                    CALL DTRSM('L','U','N','N',NR,NR,One,work(2*N+1),N,U,LDU)
//                  .. apply the permutation from the second QR factorization
                    DO 873 q = 1, NR
                        DO 872 p = 1, NR
                            work(2*N+N*NR+NR+iwork(N+p)) = U(p,q)
 872                    CONTINUE
                        DO 874 p = 1, NR
                            U(p,q) = work(2*N+N*NR+NR+p)
 874                    CONTINUE
 873                CONTINUE
                    if ( NR < N ) {
                        CALL DLASET( 'A',N-NR,NR,Zero,Zero,V(NR+1,1),LDV )
                        CALL DLASET( 'A',NR,N-NR,Zero,Zero,V(1,NR+1),LDV )
                        CALL DLASET( 'A',N-NR,N-NR,Zero,One,V(NR+1,NR+1),LDV )
                    }
                    CALL DORMQR( 'L','N',N,N,NR,work(2*N+1),N,work(N+1),
     $                      V,LDV,work(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
                } else {
//                  Last line of defense.
//#:(               This is a rather pathological case: no scaled condition
//                  improvement after two pivoted QR factorizations. Other
//                  possibility is that the rank revealing QR factorization
//                  or the condition estimator has failed, or the COND_OK
//                  is set very close to One (which is unnecessary). Normally,
//                  this branch should never be executed, but in rare cases of
//                  failure of the RRQR or condition estimator, the last line of
//                  defense ensures that DGEJSV completes the task.
//                  Compute the full SVD of L3 using DGESVJ with explicit
//                  accumulation of Jacobi rotations.
                    CALL DGESVJ( 'L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U,
     $                  LDU, work(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO )
                    SCALEM  = work(2*N+N*NR+NR+1)
                    numRank = IDNINT(work(2*N+N*NR+NR+2))
                    if ( NR < N ) {
                        CALL DLASET( 'A',N-NR,NR,Zero,Zero,V(NR+1,1),LDV )
                        CALL DLASET( 'A',NR,N-NR,Zero,Zero,V(1,NR+1),LDV )
                        CALL DLASET( 'A',N-NR,N-NR,Zero,One,V(NR+1,NR+1),LDV )
                    }
                    CALL DORMQR( 'L','N',N,N,NR,work(2*N+1),N,work(N+1),
     $                      V,LDV,work(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
//
                    CALL DORMLQ( 'L', 'T', NR, NR, NR, work(2*N+1), N,
     $                      work(2*N+N*NR+1), U, LDU, work(2*N+N*NR+NR+1),
     $                      LWORK-2*N-N*NR-NR, IERR )
                    DO 773 q = 1, NR
                        DO 772 p = 1, NR
                            work(2*N+N*NR+NR+iwork(N+p)) = U(p,q)
 772                    CONTINUE
                        DO 774 p = 1, NR
                            U(p,q) = work(2*N+N*NR+NR+p)
 774                    CONTINUE
 773                CONTINUE
//
                }
//
//              Permute the rows of V using the (column) permutation from the
//              first QRF. Also, scale the columns to make them unit in
//              Euclidean norm. This applies to all cases.
//
                TEMP1 = DSQRT(DBLE(N)) * EPSLN
                DO 1972 q = 1, N
                    DO 972 p = 1, N
                        work(2*N+N*NR+NR+iwork(p)) = V(p,q)
  972               CONTINUE
                    DO 973 p = 1, N
                        V(p,q) = work(2*N+N*NR+NR+p)
  973               CONTINUE
                    XSC = One / DNRM2( N, V(1,q), 1 )
                    if ( (XSC < (One-TEMP1)) || (XSC>(One+TEMP1)) )
     $                  CALL DSCAL( N, XSC, V(1,q), 1 )
 1972           CONTINUE
//              At this moment, V contains the right singular vectors of A.
//              Next, assemble the left singular vector matrix U (M x N).
                if ( NR < M ) {
                    CALL DLASET( 'A', M-NR, NR, Zero, Zero, U(NR+1,1), LDU )
                    if ( NR < N1 ) {
                        CALL DLASET('A',NR,N1-NR,Zero,Zero,U(1,NR+1),LDU)
                        CALL DLASET('A',M-NR,N1-NR,Zero,One,U(NR+1,NR+1),LDU)
                    }
                }
//
//              The Q matrix from the first QRF is built into the left singular
//              matrix U. This applies to all cases.
//
                CALL DORMQR( 'Left', 'No_Tr', M, N1, N, A, LDA, work, U,
     $           LDU, work(N+1), LWORK-N, IERR )

//              The columns of U are normalized. The cost is O(M*N) flops.
                TEMP1 = DSQRT(DBLE(M)) * EPSLN
                DO 1973 p = 1, NR
                    XSC = One / DNRM2( M, U(1,p), 1 )
                    if ( (XSC < (One-TEMP1)) || (XSC>(One+TEMP1)) )
     $                  CALL DSCAL( M, XSC, U(1,p), 1 )
 1973           CONTINUE
//
//              If the initial QRF is computed with row pivoting, the left
//              singular vectors must be adjusted.
//
                if ( ROWPIV )
     $              CALL DLASWP( N1, U, LDU, 1, M-1, iwork(2*N+1), -1 )
//
            } else {
//
//              .. the initial matrix A has almost orthogonal columns and
//              the second QRF is not needed
//
                CALL DLACPY( 'Upper', N, N, A, LDA, work(N+1), N )
                if ( L2PERT ) {
                    XSC = DSQRT(SMALL)
                    DO 5970 p = 2, N
                        TEMP1 = XSC * work( N + (p-1)*N + p )
                        DO 5971 q = 1, p - 1
                            work(N+(q-1)*N+p)=-DSIGN(TEMP1,work(N+(p-1)*N+q))
 5971                   CONTINUE
 5970               CONTINUE
                } else {
                    CALL DLASET( 'Lower',N-1,N-1,Zero,Zero,work(N+2),N )
                }
//
                CALL DGESVJ( 'Upper', 'U', 'N', N, N, work(N+1), N, SVA,
     $           N, U, LDU, work(N+N*N+1), LWORK-N-N*N, INFO )
//
                SCALEM  = work(N+N*N+1)
                numRank = IDNINT(work(N+N*N+2))
                DO 6970 p = 1, N
                    CALL DCOPY( N, work(N+(p-1)*N+1), 1, U(1,p), 1 )
                    CALL DSCAL( N, SVA(p), work(N+(p-1)*N+1), 1 )
 6970           CONTINUE
//
                CALL DTRSM( 'Left', 'Upper', 'NoTrans', 'No UD', N, N,
     $               One, A, LDA, work(N+1), N )
                DO 6972 p = 1, N
                    CALL DCOPY( N, work(N+p), N, V(iwork(p),1), LDV )
 6972           CONTINUE
                TEMP1 = DSQRT(DBLE(N))*EPSLN
                DO 6971 p = 1, N
                    XSC = One / DNRM2( N, V(1,p), 1 )
                    if ( (XSC < (One-TEMP1)) || (XSC>(One+TEMP1)) )
     $                  CALL DSCAL( N, XSC, V(1,p), 1 )
 6971           CONTINUE
//
//              Assemble the left singular vector matrix U (M x N).
//
                if ( N < M ) {
                    CALL DLASET( 'A',  M-N, N, Zero, Zero, U(N+1,1), LDU )
                    if ( N < N1 ) {
                        CALL DLASET( 'A',N,  N1-N, Zero, Zero,  U(1,N+1),LDU )
                        CALL DLASET( 'A',M-N,N1-N, Zero, One,U(N+1,N+1),LDU )
                    }
                }
                CALL DORMQR( 'Left', 'No Tr', M, N1, N, A, LDA, work, U,
     $           LDU, work(N+1), LWORK-N, IERR )
                TEMP1 = DSQRT(DBLE(M))*EPSLN
                DO 6973 p = 1, N1
                    XSC = One / DNRM2( M, U(1,p), 1 )
                    if ( (XSC < (One-TEMP1)) || (XSC>(One+TEMP1)) )
     $                  CALL DSCAL( M, XSC, U(1,p), 1 )
 6973           CONTINUE
//
                if ( ROWPIV )
     $              CALL DLASWP( N1, U, LDU, 1, M-1, iwork(2*N+1), -1 )
//
            }
//
//          end of the  >> almost orthogonal case <<  in the full SVD
//
        } else {
//
//          This branch deploys a preconditioned Jacobi SVD with explicitly
//          accumulated rotations. It is included as optional, mainly for
//          experimental purposes. It does perfom well, and can also be used.
//          In this implementation, this branch will be automatically activated
//          if the  condition number sigma_max(A) / sigma_min(A) is predicted
//          to be greater than the overflow threshold. This is because the
//          a posteriori computation of the singular vectors assumes robust
//          implementation of BLAS and some LAPACK procedures, capable of
//          working in presence of extreme values. Since that is not always
//          the case, ...
//
            DO 7968 p = 1, NR
                CALL DCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
 7968       CONTINUE
//
            if ( L2PERT ) {
                XSC = DSQRT(SMALL/EPSLN)
                for (IndexType q=1; q<=nr; ++q) {
                    tmp = xsc*abs(V(q,q));
                    for (IndexType p=1; p<=n; ++p) {
                        if (p>q && abs(V(p,q))<=tmp || p<q) {
                            V(p,q) = sign(tmp, V(p,q));
                        }
                        if (p<q) {
                            V(p,q) = - V(p,q);
                        }
                    }
                }
            } else {
                CALL DLASET( 'U', NR-1, NR-1, Zero, Zero, V(1,2), LDV )
            }

            CALL DGEQRF( N, NR, V, LDV, work(N+1), work(2*N+1),
     $        LWORK-2*N, IERR )
            CALL DLACPY( 'L', N, NR, V, LDV, work(2*N+1), N )
//
            DO 7969 p = 1, NR
                CALL DCOPY( NR-p+1, V(p,p), LDV, U(p,p), 1 )
 7969       CONTINUE

            if ( L2PERT ) {
                XSC = DSQRT(SMALL/EPSLN)
                DO 9970 q = 2, NR
                    DO 9971 p = 1, q - 1
                        TEMP1 = XSC * DMIN1(DABS(U(p,p)),DABS(U(q,q)))
                        U(p,q) = - DSIGN( TEMP1, U(q,p) )
 9971               CONTINUE
 9970           CONTINUE
            } else {
                CALL DLASET('U', NR-1, NR-1, Zero, Zero, U(1,2), LDU )
            }

            CALL DGESVJ( 'G', 'U', 'V', NR, NR, U, LDU, SVA,
     $        N, V, LDV, work(2*N+N*NR+1), LWORK-2*N-N*NR, INFO )
            SCALEM  = work(2*N+N*NR+1)
            numRank = IDNINT(work(2*N+N*NR+2))

            if ( NR < N ) {
                CALL DLASET( 'A',N-NR,NR,Zero,Zero,V(NR+1,1),LDV )
                CALL DLASET( 'A',NR,N-NR,Zero,Zero,V(1,NR+1),LDV )
                CALL DLASET( 'A',N-NR,N-NR,Zero,One,V(NR+1,NR+1),LDV )
            }

            CALL DORMQR( 'L','N',N,N,NR,work(2*N+1),N,work(N+1),
     $        V,LDV,work(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
//
//          Permute the rows of V using the (column) permutation from the
//          first QRF. Also, scale the columns to make them unit in
//          Euclidean norm. This applies to all cases.
//
            TEMP1 = DSQRT(DBLE(N)) * EPSLN
            DO 7972 q = 1, N
                DO 8972 p = 1, N
                    work(2*N+N*NR+NR+iwork(p)) = V(p,q)
 8972           CONTINUE
                DO 8973 p = 1, N
                    V(p,q) = work(2*N+N*NR+NR+p)
 8973           CONTINUE
                XSC = One / DNRM2( N, V(1,q), 1 )
                if ( (XSC < (One-TEMP1)) || (XSC>(One+TEMP1)) )
     $              CALL DSCAL( N, XSC, V(1,q), 1 )
 7972       CONTINUE
//
//          At this moment, V contains the right singular vectors of A.
//          Next, assemble the left singular vector matrix U (M x N).
//
            if (NR<M) {
                CALL DLASET( 'A',  M-NR, NR, Zero, Zero, U(NR+1,1), LDU )
                if (NR<N1) {
                    CALL DLASET( 'A',NR,  N1-NR, Zero, Zero,  U(1,NR+1),LDU )
                    CALL DLASET( 'A',M-NR,N1-NR, Zero, One,U(NR+1,NR+1),LDU )
                }
            }

            CALL DORMQR( 'Left', 'No Tr', M, N1, N, A, LDA, work, U,
     $        LDU, work(N+1), LWORK-N, IERR )

            if (ROWPIV)
     $         CALL DLASWP( N1, U, LDU, 1, M-1, iwork(2*N+1), -1 )


        }
        if ( TRANSP ) {
//          .. swap U and V because the procedure worked on A^t
            for (IndexType p=1; p<=n; ++p) {
                blas::swap(U(_,p),V(_,p));
            }
        }

    }
//  end of the full SVD
//
//  Undo scaling, if necessary (and possible)
//
    if ( USCAL2 <= (BIG/SVA(1))*USCAL1 ) {
        CALL DLASCL( 'G', 0, 0, USCAL1, USCAL2, NR, 1, SVA, N, IERR )
        USCAL1 = One
        USCAL2 = One
    }

    if ( NR < N ) {
        sva(_(nr+1,n)) = Zero;
    }

    work(1) = USCAL2 * SCALEM
    work(2) = USCAL1
    if ( errest ) work(3) = SCONDA
    if ( lsvec && rsvec ) {
        work(4) = condr1
        work(5) = condr2
    }
    if ( L2TRAN ) {
        work(6) = entra
        work(7) = entrat
    }

    iwork(1) = NR
    iwork(2) = numRank
    iwork(3) = warning
}
*/
    //== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VSVA, typename MU, typename MV,
          typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
jsv(JSV::Accuracy             accuracy,
    JSV::JobU                 jobU,
    JSV::JobV                 jobV,
    bool                      restrictedRange,
    bool                      considerTransA,
    bool                      perturb,
    GeMatrix<MA>              &A,
    DenseVector<VSVA>         &sva,
    GeMatrix<MU>              &U,
    GeMatrix<MV>              &V,
    DenseVector<VWORK>        &work,
    DenseVector<VIWORK>       &iwork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType  info;
    info = cxxlapack::gejsv<IndexType>(getF77Char(accuracy),
                                       getF77Char(jobU),
                                       getF77Char(jobV),
                                       restrictedRange ? 'R' : 'N',
                                       considerTransA ? 'T' : 'N',
                                       perturb ? 'P' : 'N',
                                       A.numRows(),
                                       A.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       sva.data(),
                                       U.data(),
                                       U.leadingDimension(),
                                       V.data(),
                                       V.leadingDimension(),
                                       work.data(),
                                       work.length(),
                                       iwork.data());
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename VSVA, typename MU, typename MV,
          typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
jsv(JSV::Accuracy             accuracy,
    JSV::JobU                 jobU,
    JSV::JobV                 jobV,
    bool                      restrictedRange,
    bool                      considerTransA,
    bool                      perturb,
    GeMatrix<MA>              &A,
    DenseVector<VSVA>         &sva,
    GeMatrix<MU>              &U,
    GeMatrix<MV>              &V,
    DenseVector<VWORK>        &work,
    DenseVector<VIWORK>       &iwork)
{
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

    ASSERT(sva.firstIndex()==1);
    ASSERT(sva.length()==n);

    ASSERT(U.firstRow()==1);
    ASSERT(U.firstCol()==1);

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);

    ASSERT(iwork.length()==m+3*n);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView        A_org     = A;
    typename DenseVector<VSVA>::NoView   sva_org   = sva;
    typename GeMatrix<MU>::NoView        U_org     = U;
    typename GeMatrix<MV>::NoView        V_org     = V;
    typename DenseVector<VWORK>::NoView  work_org  = work;
#   endif

//
//  Call implementation
//
    /*
    IndexType info = jsv_generic(accuracy, jobU, jobV,
                                 restrictedRange, considerTransA, perturb,
                                 A, sva, U, V, work);
    */
    IndexType info = external::jsv(accuracy, jobU, jobV,
                                   restrictedRange, considerTransA, perturb,
                                   A, sva, U, V, work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView        A_generic     = A;
    typename DenseVector<VSVA>::NoView   sva_generic   = sva;
    typename GeMatrix<MU>::NoView        U_generic     = U;
    typename GeMatrix<MV>::NoView        V_generic     = V;
    typename DenseVector<VWORK>::NoView  work_generic  = work;
//
//  restore output arguments
//
    A    = A_org;
    sva  = sva_org;
    U    = U_org;
    V    = V_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::jsv(accuracy, jobU, jobV,
                                    restrictedRange, considerTransA, perturb,
                                    A, sva, U, V, work, iwork);
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
    if (! isIdentical(U_generic, U, "U_generic", "U")) {
        std::cerr << "CXXLAPACK: U_generic = " << U_generic << std::endl;
        std::cerr << "F77LAPACK: U = " << U << std::endl;
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
        std::cerr << "error in: jsv.tcc " << std::endl;
        ASSERT(0);
    } else {
        std::cerr << "passed: jsv.tcc " << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VSVA, typename MU, typename MV,
         typename VWORK, typename VIWORK>
typename MA::IndexType
jsv(JSV::Accuracy             accuracy,
    JSV::JobU                 jobU,
    JSV::JobV                 jobV,
    bool                      restrictedRange,
    bool                      considerTransA,
    bool                      perturb,
    MA                        &&A,
    VSVA                      &&sva,
    MU                        &&U,
    MV                        &&V,
    VWORK                     &&work,
    VIWORK                    &&iwork)
{
    typename MA::IndexType info;

    CHECKPOINT_ENTER;
    info = jsv(accuracy, jobU, jobV,
               restrictedRange, considerTransA, perturb,
               A, sva, U, V, work, iwork);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SVD_JSV_TCC
