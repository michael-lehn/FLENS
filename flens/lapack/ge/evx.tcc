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
      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
     $                   VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
      SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     $                   LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE,
     $                   RCONDV, WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_EVX_TCC
#define FLENS_LAPACK_GE_EVX_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)evx_wsq [worksize query, real variant] --------------------------------

template <typename MA>
Pair<typename GeMatrix<MA>::IndexType>
evx_wsq_impl(bool               computeVL,
             bool               computeVR,
             SENSE::Sense       sense,
             const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType n = A.numRows();
    const bool wantSN = (sense==SENSE::None);
    const bool wantSE = (sense==SENSE::EigenvaluesOnly);

    IndexType minWork, maxWork;

    if (n==0) {
        minWork = 1;
        maxWork = 1;
    } else {
        maxWork = n + n*ilaenv<T>(1, "GEHRD", "", n, 1, n, 0);

        IndexType hseqrWork;
        if (computeVL) {
            hseqrWork = hseqr_wsq(HSEQR::Schur, HSEQR::NoInit,
                                  IndexType(1), n, A);
        } else if (computeVL) {
            hseqrWork = hseqr_wsq(HSEQR::Schur, HSEQR::NoInit,
                                  IndexType(1), n, A);
        } else {
            if (wantSN) {
                hseqrWork = hseqr_wsq(HSEQR::Eigenvalues, HSEQR::No,
                                      IndexType(1), n, A);
            } else {
                hseqrWork = hseqr_wsq(HSEQR::Schur, HSEQR::No,
                                      IndexType(1), n, A);
            }
        }

        if ((!computeVL) && (!computeVR)) {
            minWork = 2*n;
            if (!wantSN) {
                minWork = max(minWork, n*n+6*n);
            }
            maxWork = max(maxWork, hseqrWork);
            if (!wantSN) {
                maxWork = max(maxWork, n*n+6*n);
            }
        } else {
            minWork = 3*n;
            if ((!wantSN) && (!wantSE)) {
                minWork = max(minWork, n*n + 6*n);
            }
            maxWork = max(maxWork, hseqrWork);
            maxWork = max(maxWork, n+(n-1*ilaenv<T>(1, "ORGHR", "", n, 1, n)));
            if ((!wantSN) && (!wantSE)) {
                maxWork = max(maxWork, n*n + 6*n);
            }
            maxWork = max(maxWork, 3*n);
        }
        maxWork = max(maxWork, minWork);
    }

    return Pair<typename GeMatrix<MA>::IndexType>(minWork, maxWork);
}

//-- (ge)evx [real variant] ----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename RCONDE, typename RCONDV, typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
evx_impl(BALANCE::Balance     balance,
         bool                 computeVL,
         bool                 computeVR,
         SENSE::Sense         sense,
         GeMatrix<MA>         &A,
         DenseVector<VWR>     &wr,
         DenseVector<VWI>     &wi,
         GeMatrix<MVL>        &VL,
         GeMatrix<MVR>        &VR,
         IndexType            &iLo,
         IndexType            &iHi,
         DenseVector<VSCALE>  &scale,
         ABNORM               &ABNorm,
         DenseVector<RCONDE>  &rCondE,
         DenseVector<RCONDV>  &rCondV,
         DenseVector<VWORK>   &work,
         DenseVector<VIWORK>  &iWork)
{
    using cxxblas::pow;

    typedef typename GeMatrix<MA>::ElementType   T;

    const IndexType             n = A.numRows();
    const Underscore<IndexType> _;

    const T   Zero(0), One(1);
//
//  .. Local Arrays ..
//
    bool selectData_[1];
    DenseVectorView<bool>
        select = typename DenseVectorView<bool>::Engine(1, selectData_);
//
//  Test the input arguments
//
    IndexType info = 0;

    const bool wantSN = (sense==SENSE::None);
    const bool wantSV = (sense==SENSE::InvariantSubspaceOnly);
    const bool wantSB = (sense==SENSE::Both);
//
//  Compute workspace
//   (Note: Comments in the code beginning "Workspace:" describe the
//    minimal amount of workspace needed at that point in the code,
//    as well as the preferred amount for good performance.
//    NB refers to the optimal block size for the immediately
//    following subroutine, as returned by ILAENV.
//    HSWORK refers to the workspace preferred by DHSEQR, as
//    calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
//    the worst case.)
//
    Pair<IndexType> wsQuery = evx_wsq(computeVL, computeVR, sense, A);
    IndexType minWork = wsQuery.first;
    IndexType maxWork = wsQuery.second;

    if (work.length()!=0 && work.length()<minWork) {
        ASSERT(0);
    } else if (work.length()==0) {
        work.resize(maxWork);
    }
    work(1) = maxWork;

//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Get machine constants
//
    const T eps = lamch<T>(Precision);
    T smallNum = lamch<T>(SafeMin);
    T bigNum = One / smallNum;
    labad(smallNum, bigNum);
    smallNum = sqrt(smallNum) / eps;
    bigNum = One / smallNum;
//
//  Scale A if max element outside range [SMLNUM,BIGNUM]
//
    IndexType iCond = 0;
    const T ANorm = lan(MaximumNorm, A);
    bool scaleA = false;
    T cScale;
    if (ANorm>Zero && ANorm<smallNum) {
        scaleA = true;
        cScale = smallNum;
    } else if (ANorm>bigNum) {
        scaleA = true;
        cScale = bigNum;
    }
    if (scaleA) {
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0), ANorm, cScale, A);
    }
//
//  Balance the matrix and compute ABNRM
//
    bal(balance, A, iLo, iHi, scale);
    ABNorm = lan(OneNorm, A);
    if (scaleA) {
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
              cScale, ANorm, ABNorm);
    }
//
//  Reduce to upper Hessenberg form
//  (Workspace: need 2*N, prefer N+N*NB)
//
    IndexType iTau = 1;
    IndexType iWrk = iTau + n;
    IndexType lWork = work.length();

    auto tau = work(_(iTau, iTau+n-2));
    auto hrdWork = work(_(iWrk, lWork));

    hrd(iLo, iHi, A, tau, hrdWork);

    if (computeVL) {
//
//      Want left eigenvectors
//      Copy Householder vectors to VL
//
        VL.lower() = A.lower();
//
//      Generate orthogonal matrix in VL
//      (Workspace: need 2*N-1, prefer N+(N-1)*NB)
//
        orghr(iLo, iHi, VL, tau, hrdWork);
//
//      Perform QR iteration, accumulating Schur vectors in VL
//      (Workspace: need 1, prefer HSWORK (see comments) )
//
        iWrk = iTau;
        auto hseqrWork = work(_(iWrk, lWork));
        info = hseqr(HSEQR::Schur, HSEQR::NoInit, iLo, iHi, A,
                     wr, wi, VL, hseqrWork);

        if (computeVR) {
//
//          Want left and right eigenvectors
//          Copy Schur vectors to VR
//
            VR = VL;
        }

    } else if (computeVR) {
//
//      Want right eigenvectors
//      Copy Householder vectors to VR
//
        VR.lower() = A.lower();
//
//      Generate orthogonal matrix in VR
//      (Workspace: need 2*N-1, prefer N+(N-1)*NB)
//
        orghr(iLo, iHi, VR, tau, hrdWork);
//
//      Perform QR iteration, accumulating Schur vectors in VR
//      (Workspace: need 1, prefer HSWORK (see comments) )
//
        iWrk = iTau;
        auto hseqrWork = work(_(iWrk, lWork));
        info = hseqr(HSEQR::Schur, HSEQR::NoInit, iLo, iHi, A,
                     wr, wi, VR, hseqrWork);

    } else {
//
//      Compute eigenvalues only
//      If condition numbers desired, compute Schur form
//
        HSEQR::Job  job = (wantSN) ? HSEQR::Eigenvalues
                                   : HSEQR::Schur;
//
//      (Workspace: need 1, prefer HSWORK (see comments) )
//
        iWrk = iTau;
        auto hseqrWork = work(_(iWrk, lWork));
        info = hseqr(job, HSEQR::No, iLo, iHi, A, wr, wi, VR, hseqrWork);
    }
//
//  If INFO > 0 from DHSEQR, then quit
//
    if (info==0) {

        if (computeVL || computeVR) {
//
//          Compute left and/or right eigenvectors
//          (Workspace: need 3*N)
//
            IndexType nOut;
            auto      trevcWork = work(_(iWrk, lWork));

            trevc(computeVL, computeVR, TREVC::Backtransform, select,
                  A, VL, VR, n, nOut, trevcWork);
        }
//
//      Compute condition numbers if desired
//      (Workspace: need N*N+6*N unless SENSE = 'E')
//
        if (!wantSN) {
            IndexType nOut;

            IndexType n_ = (sense!=SENSE::EigenvaluesOnly) ? n : 0;
            GeMatrixView<T>  Work(n_, n+6, work(_(iWrk, lWork)), n);

            trsna(TRSNA::Job(sense), TRSNA::All, select, A, VL, VR,
                  rCondE, rCondV, n, nOut, Work, iWork);
        }

        if (computeVL) {
//
//          Undo balancing of left eigenvectors
//
            bak(balance, Left, iLo, iHi, scale, VL);
//
//          Normalize left eigenvectors and make largest component real
//
            for (IndexType i=1; i<=n; ++i) {
                if (wi(i)==Zero) {
                    VL(_,i) *= One / blas::nrm2(VL(_,i));
                } else if (wi(i)>Zero) {
                    const T scl = One / lapy2(blas::nrm2(VL(_,i)),
                                              blas::nrm2(VL(_,i+1)));
                    VL(_,i)     *= scl;
                    VL(_,i+1)   *= scl;
                    for (IndexType k=1; k<=n; ++k) {
                        work(k) = pow(VL(k,i),2) + pow(VL(k,i+1),2);
                    }
                    IndexType k = blas::iamax(work(_(1,n)));
                    T cs, sn, r;
                    lartg(VL(k,i), VL(k,i+1), cs, sn, r);
                    blas::rot(VL(_,i), VL(_,i+1), cs, sn);
                    VL(k,i+1) = Zero;
                }
            }
        }

        if (computeVR) {
//
//          Undo balancing of right eigenvectors
//
            bak(BALANCE::Both, Right, iLo, iHi, scale, VR);
//
//          Normalize right eigenvectors and make largest component real
//
            for (IndexType i=1; i<=n; ++i) {
                if (wi(i)==Zero) {
                    VR(_,i) *= One / blas::nrm2(VR(_,i));
                } else if (wi(i)>Zero) {
                    const T scl = One / lapy2(blas::nrm2(VR(_,i)),
                                              blas::nrm2(VR(_,i+1)));
                    VR(_,i)     *= scl;
                    VR(_,i+1)   *= scl;
                    for (IndexType k=1; k<=n; ++k) {
                        work(k) = pow(VR(k,i), 2) + pow(VR(k,i+1), 2);
                    }
                    IndexType k = blas::iamax(work(_(1,n)));
                    T cs, sn, r;
                    lartg(VR(k,i), VR(k,i+1), cs, sn, r);
                    blas::rot(VR(_,i), VR(_,i+1), cs, sn);
                    VR(k,i+1) = Zero;
                }
            }
        }
    }

//
//  Undo scaling if necessary
//
    if (scaleA) {
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
              cScale, ANorm, wr(_(info+1,n)));
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
              cScale, ANorm, wi(_(info+1,n)));

        if (info==0) {
            if ((wantSV || wantSB) && iCond==0) {
                lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                      cScale, ANorm, rCondV);
            }
        } else {
            lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                  cScale, ANorm, wr(_(1,iLo-1)));
            lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                  cScale, ANorm, wi(_(1,iLo-1)));
        }
    }

    work(1) = maxWork;
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)evx_wsq [worksize query, real variant] --------------------------------

template <typename MA>
typename RestrictTo<IsNotComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
evx_wsq_impl(bool               computeVL,
             bool               computeVR,
             SENSE::Sense       sense,
             const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const bool wantSN = (sense==SENSE::None);
    const bool wantSE = (sense==SENSE::EigenvaluesOnly);

//
//  Compute minimal worksize
//
    IndexType  n = A.numRows();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        if ((!computeVL) && (!computeVR)) {
            minWork = 2*n;
            if (!wantSN) {
                minWork = max(minWork, n*n+6*n);
            }
        } else {
            minWork = 3*n;
            if ((!wantSN) && (!wantSE)) {
                minWork = max(minWork, n*n + 6*n);
            }
        }
    }

//
//  Compute optimal worksize
//
    const IndexType     LDVL   = A.numRows();
    const IndexType     LDVR   = A.numRows();
    IndexType           IDUMMY;
    T                   DUMMY;
    const IndexType     LWORK  = -1;
    T                   WORK;

    cxxlapack::geevx<IndexType>('N',
                                computeVL ? 'V' : 'N',
                                computeVR ? 'V' : 'N',
                                getF77Char(sense),
                                A.numRows(),
                                &DUMMY,
                                A.leadingDimension(),
                                &DUMMY,
                                &DUMMY,
                                &DUMMY,
                                LDVL,
                                &DUMMY,
                                LDVR,
                                IDUMMY,
                                IDUMMY,
                                &DUMMY,
                                DUMMY,
                                &DUMMY,
                                &DUMMY,
                                &WORK,
                                LWORK,
                                &IDUMMY);
    return Pair<IndexType>(minWork, WORK);
}

//-- (ge)evx_wsq [worksize query, complex variant] -----------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
evx_wsq_impl(bool               computeVL,
             bool               computeVR,
             SENSE::Sense       sense,
             const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     RT;
    typedef typename GeMatrix<MA>::IndexType            IndexType;

    const bool wantSN = (sense==SENSE::None);
    const bool wantSE = (sense==SENSE::EigenvaluesOnly);

//
//  Compute minimal worksize
//
    IndexType  n = A.numRows();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        minWork = 2*n;
        if (!wantSN && !wantSE) {
            minWork = max(minWork, n*n+2*n);
        }
    }

//
//  Compute optimal worksize
//
    const IndexType     LDVL   = A.numRows();
    const IndexType     LDVR   = A.numRows();
    IndexType           IDUMMY;
    T                   DUMMY;
    RT                  RDUMMY;
    const IndexType     LWORK  = -1;
    T                   WORK;

    cxxlapack::geevx<IndexType>('N',
                                computeVL ? 'V' : 'N',
                                computeVR ? 'V' : 'N',
                                getF77Char(sense),
                                A.numRows(),
                                &DUMMY,
                                A.leadingDimension(),
                                &DUMMY,                     // w
                                &DUMMY,                     // VL
                                LDVL,                       // ldVL
                                &DUMMY,                     // VR
                                LDVR,                       // ldVR
                                IDUMMY,                     // iLo
                                IDUMMY,                     // iHi
                                &RDUMMY,                    // scale
                                RDUMMY,                     // ABnorm
                                &RDUMMY,                    // rCondE
                                &RDUMMY,                    // rCondV
                                &WORK,                      // work
                                LWORK,                      // lWork
                                &RDUMMY);                   // rWork
    return Pair<IndexType>(minWork, WORK.real());
}

//-- (ge)evx [real variant] ----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename RCONDE, typename RCONDV, typename VWORK, typename VIWORK>
typename GeMatrix<MA>::IndexType
evx_impl(BALANCE::Balance     balance,
         bool                 computeVL,
         bool                 computeVR,
         SENSE::Sense         sense,
         GeMatrix<MA>         &A,
         DenseVector<VWR>     &wr,
         DenseVector<VWI>     &wi,
         GeMatrix<MVL>        &VL,
         GeMatrix<MVR>        &VR,
         IndexType            &iLo,
         IndexType            &iHi,
         DenseVector<VSCALE>  &scale,
         ABNORM               &abNorm,
         DenseVector<RCONDE>  &rCondE,
         DenseVector<RCONDV>  &rCondV,
         DenseVector<VWORK>   &work,
         DenseVector<VIWORK>  &iWork)
{
    if (work.length()==0) {
        const auto ws = evx_wsq_impl(computeVL, computeVR, sense, A);
        work.resize(ws.second, 1);
    }

    IndexType  info;
    info = cxxlapack::geevx<IndexType>(getF77Char(balance),
                                       computeVL ? 'V' : 'N',
                                       computeVR ? 'V' : 'N',
                                       getF77Char(sense),
                                       A.numRows(),
                                       A.data(),
                                       A.leadingDimension(),
                                       wr.data(),
                                       wi.data(),
                                       VL.data(),
                                       VL.leadingDimension(),
                                       VR.data(),
                                       VR.leadingDimension(),
                                       iLo,
                                       iHi,
                                       scale.data(),
                                       abNorm,
                                       rCondE.data(),
                                       rCondV.data(),
                                       work.data(),
                                       work.length(),
                                       iWork.data());
    ASSERT(info>=0);
    return info;
}

//-- (ge)evx [complex variant] -------------------------------------------------

template <typename MA, typename VW, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename RCONDE, typename RCONDV, typename VWORK, typename VRWORK>
typename GeMatrix<MA>::IndexType
evx_impl(BALANCE::Balance     balance,
         bool                 computeVL,
         bool                 computeVR,
         SENSE::Sense         sense,
         GeMatrix<MA>         &A,
         DenseVector<VW>      &w,
         GeMatrix<MVL>        &VL,
         GeMatrix<MVR>        &VR,
         IndexType            &iLo,
         IndexType            &iHi,
         DenseVector<VSCALE>  &scale,
         ABNORM               &abNorm,
         DenseVector<RCONDE>  &rCondE,
         DenseVector<RCONDV>  &rCondV,
         DenseVector<VWORK>   &work,
         DenseVector<VRWORK>  &rWork)
{
    if (work.length()==0) {
        const auto ws = evx_wsq_impl(computeVL, computeVR, sense, A);
        work.resize(ws.second, 1);
    }

    IndexType  info;
    info = cxxlapack::geevx<IndexType>(getF77Char(balance),
                                       computeVL ? 'V' : 'N',
                                       computeVR ? 'V' : 'N',
                                       getF77Char(sense),
                                       A.numRows(),
                                       A.data(),
                                       A.leadingDimension(),
                                       w.data(),
                                       VL.data(),
                                       VL.leadingDimension(),
                                       VR.data(),
                                       VR.leadingDimension(),
                                       iLo,
                                       iHi,
                                       scale.data(),
                                       abNorm,
                                       rCondE.data(),
                                       rCondV.data(),
                                       work.data(),
                                       work.length(),
                                       rWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)evx [real variant] ----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename VRCONDE, typename VRCONDV, typename VWORK, typename VIWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MVL>::value
                 && IsRealGeMatrix<MVR>::value
                 && IsInteger<IndexType>::value
                 && IsRealDenseVector<VSCALE>::value
                 && IsNotComplex<ABNORM>::value
                 && IsRealDenseVector<VRCONDE>::value
                 && IsRealDenseVector<VRCONDV>::value
                 && IsRealDenseVector<VWORK>::value
                 && IsIntegerDenseVector<VIWORK>::value,
         IndexType>::Type
evx(BALANCE::Balance     balance,
    bool                 computeVL,
    bool                 computeVR,
    SENSE::Sense         sense,
    MA                   &&A,
    VWR                  &&wr,
    VWI                  &&wi,
    MVL                  &&VL,
    MVR                  &&VR,
    IndexType            &iLo,
    IndexType            &iHi,
    VSCALE               &&scale,
    ABNORM               &abNorm,
    VRCONDE              &&rCondE,
    VRCONDV              &&rCondV,
    VWORK                &&work,
    VIWORK               &&iWork)
{
    LAPACK_DEBUG_OUT("(ge)evx [real]");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<VWR>::Type       VectorWR;
    typedef typename RemoveRef<VWI>::Type       VectorWI;
    typedef typename RemoveRef<MVL>::Type       MatrixVL;
    typedef typename RemoveRef<MVR>::Type       MatrixVR;
    typedef typename RemoveRef<VSCALE>::Type    VectorScale;
    typedef typename RemoveRef<VRCONDE>::Type   VectorRCondE;
    typedef typename RemoveRef<VRCONDV>::Type   VectorRCondV;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
    typedef typename RemoveRef<VIWORK>::Type    VectorIWork;
#   endif

    const IndexType n = A.numRows();

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);

    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==0 || wr.length()==n);

    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==0 || wi.length()==n);

    if (computeVL) {
        ASSERT(VL.numRows()==VL.numCols());
        ASSERT(VL.numRows()==0 || VL.numRows()==n);
        ASSERT(VL.firstRow()==1);
        ASSERT(VL.firstCol()==1);
    }

    if (computeVR) {
        ASSERT(VR.numRows()==VR.numCols());
        ASSERT(VR.numRows()==0 || VR.numRows()==n);
        ASSERT(VR.firstRow()==1);
        ASSERT(VR.firstCol()==1);
    }

    ASSERT(scale.firstIndex()==1);
    ASSERT(scale.length()==0 || scale.length()==n);

    ASSERT(rCondE.firstIndex()==1);
    ASSERT(rCondE.length()==0 || rCondE.length()==n);

    ASSERT(rCondV.firstIndex()==1);
    ASSERT(rCondV.length()==0 || rCondV.length()==n);

    ASSERT(iWork.length()==0 || iWork.length()==2*(n-1));
#   endif

//
//  Resize output arguments if they are empty and needed
//
    if (wr.length()==0) {
        wr.resize(n, 1);
    }
    if (wi.length()==0) {
        wi.resize(n, 1);
    }
    if (computeVL && VL.numRows()==0) {
        VL.resize(n, n, 1, 1);
    }
    if (computeVR && VR.numRows()==0) {
        VR.resize(n, n, 1, 1);
    }
    if (scale.length()==0) {
        scale.resize(n, 1);
    }
    if (rCondE.length()==0) {
        rCondE.resize(n, 1);
    }
    if (rCondV.length()==0) {
        rCondV.resize(n, 1);
    }
    if (iWork.length()==0) {
        iWork.resize(2*(n-1), 1);
    }
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_org       = A;
    typename VectorWR::NoView       wr_org      = wr;
    typename VectorWI::NoView       wi_org      = wi;
    typename MatrixVL::NoView       VL_org      = VL;
    typename MatrixVR::NoView       VR_org      = VR;
    IndexType                       iLo_org     = iLo;
    IndexType                       iHi_org     = iHi;
    typename VectorScale::NoView    scale_org   = scale;
    ABNORM                          abNorm_org  = abNorm;
    typename VectorRCondE::NoView   rCondE_org  = rCondE;
    typename VectorRCondV::NoView   rCondV_org  = rCondV;
    typename VectorWork::NoView     work_org    = work;
    typename VectorIWork::NoView    iWork_org   = iWork;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::evx_impl(balance, computeVL, computeVR,
                                             sense, A, wr, wi, VL, VR,
                                             iLo, iHi, scale, abNorm,
                                             rCondE, rCondV, work, iWork);
#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView        A_generic       = A;
    typename VectorWR::NoView       wr_generic      = wr;
    typename VectorWI::NoView       wi_generic      = wi;
    typename MatrixVL::NoView       VL_generic      = VL;
    typename MatrixVR::NoView       VR_generic      = VR;
    IndexType                       iLo_generic     = iLo;
    IndexType                       iHi_generic     = iHi;
    typename VectorScale::NoView    scale_generic   = scale;
    ABNORM                          abNorm_generic  = abNorm;
    typename VectorRCondE::NoView   rCondE_generic  = rCondE;
    typename VectorRCondV::NoView   rCondV_generic  = rCondV;
    typename VectorWork::NoView     work_generic    = work;
    typename VectorIWork::NoView    iWork_generic   = iWork;
//
//  restore output arguments
//
    A       = A_org;
    wr      = wr_org;
    wi      = wi_org;
    VL      = VL_org;
    VR      = VR_org;
    iLo     = iLo_org;
    iHi     = iHi_org;
    scale   = scale_org;
    abNorm  = abNorm_org;
    rCondE  = rCondE_org;
    rCondV  = rCondV_org;
    work    = work_org;
    iWork   = iWork_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType info_ = external::evx_impl(balance, computeVL, computeVR,
                                         sense, A, wr, wi, VL, VR, iLo, iHi,
                                         scale, abNorm, rCondE, rCondV,
                                         work, iWork);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(wr_generic, wr, "wr_generic", "wr")) {
        std::cerr << "CXXLAPACK: wr_generic = " << wr_generic << std::endl;
        std::cerr << "F77LAPACK: wr = " << wr << std::endl;
        failed = true;
    }

    if (! isIdentical(wi_generic, wi, "wi_generic", "wi")) {
        std::cerr << "CXXLAPACK: wi_generic = " << wi_generic << std::endl;
        std::cerr << "F77LAPACK: wi = " << wi << std::endl;
        failed = true;
    }

    if (! isIdentical(VL_generic, VL, "VL_generic", "VL")) {
        std::cerr << "CXXLAPACK: VL_generic = " << VL_generic << std::endl;
        std::cerr << "F77LAPACK: VL = " << VL << std::endl;
        failed = true;
    }

    if (! isIdentical(VR_generic, VR, "VR_generic", "VR")) {
        std::cerr << "CXXLAPACK: VR_generic = " << VR_generic << std::endl;
        std::cerr << "F77LAPACK: VR = " << VR << std::endl;
        failed = true;
    }

    if (! isIdentical(iLo_generic, iLo, "iLo_generic", "iLo")) {
        std::cerr << "CXXLAPACK: iLo_generic = " << iLo_generic << std::endl;
        std::cerr << "F77LAPACK: iLo = " << iLo << std::endl;
        failed = true;
    }

    if (! isIdentical(iHi_generic, iHi, "iHi_generic", "iHi")) {
        std::cerr << "CXXLAPACK: iHi_generic = " << iHi_generic << std::endl;
        std::cerr << "F77LAPACK: iHi = " << iHi << std::endl;
        failed = true;
    }

    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = " << scale << std::endl;
        failed = true;
    }

    if (! isIdentical(abNorm_generic, abNorm, "abNorm_generic", "abNorm")) {
        std::cerr << "CXXLAPACK: abNorm_generic = "
                  << abNorm_generic << std::endl;
        std::cerr << "F77LAPACK: abNorm = " << abNorm << std::endl;
        failed = true;
    }

    if (! isIdentical(rCondE_generic, rCondE, "rCondE_generic", "rCondE")) {
        std::cerr << "CXXLAPACK: rCondE_generic = "
                  << rCondE_generic << std::endl;
        std::cerr << "F77LAPACK: rCondE = " << rCondE << std::endl;
        failed = true;
    }

    if (! isIdentical(rCondV_generic, rCondV, "rCondV_generic", "rCondV")) {
        std::cerr << "CXXLAPACK: rCondV_generic = "
                  << rCondV_generic << std::endl;
        std::cerr << "F77LAPACK: rCondV = " << rCondV << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(iWork_generic, iWork, "iWork_generic", "iWork")) {
        std::cerr << "CXXLAPACK: iWork_generic = "
                  << iWork_generic << std::endl;
        std::cerr << "F77LAPACK: iWork = " << iWork << std::endl;
        failed = true;
    }

    if (! isIdentical(info, info_, " info", "info_")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: info_ = " << info_ << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
//        std::cerr << "passed: evx.tcc" << std::endl;
    }
#   endif

    return info;
}


#ifdef USE_CXXLAPACK

//-- (ge)evx [complex variant] -------------------------------------------------

template <typename MA, typename VW, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename VRCONDE, typename VRCONDV, typename VWORK, typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MVL>::value
                 && IsComplexGeMatrix<MVR>::value
                 && IsInteger<IndexType>::value
                 && IsRealDenseVector<VSCALE>::value
                 && IsNotComplex<ABNORM>::value
                 && IsRealDenseVector<VRCONDE>::value
                 && IsRealDenseVector<VRCONDV>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         IndexType>::Type
evx(BALANCE::Balance     balance,
    bool                 computeVL,
    bool                 computeVR,
    SENSE::Sense         sense,
    MA                   &&A,
    VW                   &&w,
    MVL                  &&VL,
    MVR                  &&VR,
    IndexType            &iLo,
    IndexType            &iHi,
    VSCALE               &&scale,
    ABNORM               &abNorm,
    VRCONDE              &&rCondE,
    VRCONDV              &&rCondV,
    VWORK                &&work,
    VRWORK               &&rWork)
{
    LAPACK_DEBUG_OUT("(ge)evx [complex]");

    const IndexType n = A.numRows();

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==0 || w.length()==n);

    if (computeVL) {
        ASSERT(VL.numRows()==VL.numCols());
        ASSERT(VL.numRows()==0 || VL.numRows()==n);
        ASSERT(VL.firstRow()==1);
        ASSERT(VL.firstCol()==1);
    }

    if (computeVR) {
        ASSERT(VR.numRows()==VR.numCols());
        ASSERT(VR.numRows()==0 || VR.numRows()==n);
        ASSERT(VR.firstRow()==1);
        ASSERT(VR.firstCol()==1);
    }

    ASSERT(scale.firstIndex()==1);
    ASSERT(scale.length()==0 || scale.length()==n);

    ASSERT(rCondE.firstIndex()==1);
    ASSERT(rCondE.length()==0 || rCondE.length()==n);

    ASSERT(rCondV.firstIndex()==1);
    ASSERT(rCondV.length()==0 || rCondV.length()==n);

    ASSERT(rWork.length()==0 || rWork.length()==2*n);
#   endif

//
//  Resize output arguments if they are empty (and needed).
//
    if (w.length()==0) {
        w.resize(n, 1);
    }
    if (computeVL && VL.numRows()==0) {
        VL.resize(n, n, 1, 1);
    }
    if (computeVR && VR.numRows()==0) {
        VR.resize(n, n, 1, 1);
    }
    if (scale.length()==0) {
        scale.resize(n, 1);
    }
    if (rCondE.length()==0) {
        rCondE.resize(n, 1);
    }
    if (rCondV.length()==0) {
        rCondV.resize(n, 1);
    }
    if (rWork.length()==0) {
        rWork.resize(2*n, 1);
    }

//
//  Call external implementation
//
    IndexType info = external::evx_impl(balance, computeVL, computeVR,
                                        sense, A, w, VL, VR,
                                        iLo, iHi, scale, abNorm,
                                        rCondE, rCondV, work, rWork);
    return info;
}

#endif // USE_CXXLAPACK

//-- (ge)evx_wsq [worksize query, real variant] --------------------------------

template <typename MA>
typename RestrictTo<IsRealGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
evx_wsq(bool computeVL, bool computeVR, SENSE::Sense sense, const MA &A)
{
    return LAPACK_SELECT::evx_wsq_impl(computeVL, computeVR, sense, A);
}

//-- (ge)evx_wsq [worksize query, complex variant] -----------------------------

#ifdef USE_CXXLAPACK

template <typename MA>
typename RestrictTo<IsComplexGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
evx_wsq(bool computeVL, bool computeVR, SENSE::Sense sense, const MA &A)
{
    return external::evx_wsq_impl(computeVL, computeVR, sense, A);
}

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_EVX_TCC
