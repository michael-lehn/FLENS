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
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_EV_TCC
#define FLENS_LAPACK_IMPL_EV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)ev_wsq [worksize query, real variant] ---------------------------------

template <typename MA>
Pair<typename GeMatrix<MA>::IndexType>
ev_wsq_impl(bool computeVL, bool computeVR, const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType n = A.numRows();

    IndexType minWork, maxWork;

    if (n==0) {
        minWork = 1;
        maxWork = 1;
    } else {
        maxWork = 2*n + n*ilaenv<T>(1, "GEHRD", "", n, 1, n, 0);
        if (computeVL) {
            minWork = 4*n;
            maxWork = max(maxWork,
                          2*n + (n-1)*ilaenv<T>(1, "ORGHR", "", n, 1, n));
            IndexType hseqrWork = hseqr_wsq(HSEQR::Schur, HSEQR::NoInit,
                                            IndexType(1), n, A);
            maxWork = max(max(maxWork, n+1), n+hseqrWork);
            maxWork = max(maxWork, 4*n);
        } else if (computeVR) {
            minWork = 4*n;
            maxWork = max(maxWork,
                          2*n + (n-1)*ilaenv<T>(1, "ORGHR", "", n, 1, n));
            IndexType hseqrWork = hseqr_wsq(HSEQR::Schur, HSEQR::NoInit,
                                            IndexType(1), n, A);
            maxWork = max(max(maxWork, n+1), n+hseqrWork);
            maxWork = max(maxWork, 4*n);
        } else {
            minWork = 3*n;
            IndexType hseqrWork = hseqr_wsq(HSEQR::Eigenvalues, HSEQR::No,
                                            IndexType(1), n, A);
            maxWork = max(max(maxWork, n+1), n+hseqrWork);
        }
        maxWork = max(maxWork, minWork);
    }
    return Pair<IndexType>(minWork, maxWork);
}

//-- (ge)ev [real variant] -----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
typename GeMatrix<MA>::IndexType
ev_impl(bool                computeVL,
        bool                computeVR,
        GeMatrix<MA>        &A,
        DenseVector<VWR>    &wr,
        DenseVector<VWI>    &wi,
        GeMatrix<MVL>       &VL,
        GeMatrix<MVR>       &VR,
        DenseVector<VWORK>  &work)
{
    using flens::pow;
    using std::sqrt;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const T Zero(0), One(1);
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
    Pair<IndexType> wsQuery = ev_wsq(computeVL, computeVR, A);
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
        return 0;
    }
//
//  Get machine constants
//
    const T Eps = lamch<T>(Precision);
    T smallNum = lamch<T>(SafeMin);
    T bigNum = One / smallNum;
    labad(smallNum, bigNum);
    smallNum = sqrt(smallNum) / Eps;
    bigNum = One / smallNum;
//
//  Scale A if max element outside range [SMLNUM,BIGNUM]
//
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
        lascl(LASCL::FullMatrix, 0, 0, ANorm, cScale, A);
    }
//
//  Balance the matrix
//  (Workspace: need N)
//
    IndexType iBal = 1;
    IndexType iLo, iHi;

    auto balWork = work(_(iBal,iBal+n-1));
    bal(BALANCE::Both, A, iLo, iHi, balWork);

//
//  Reduce to upper Hessenberg form
//  (Workspace: need 3*N, prefer 2*N+N*NB)
//
    IndexType iTau = iBal + n;
    IndexType iWork = iTau + n;
    IndexType lWork = work.length();

    auto tau = work(_(iTau, iTau+n-2));
    auto hrdWork = work(_(iWork, lWork));

    hrd(iLo, iHi, A, tau, hrdWork);

    IndexType info = 0;
    if (computeVL) {
//
//      Want left eigenvectors
//      Copy Householder vectors to VL
//
        VL.lower() = A.lower();
//
//      Generate orthogonal matrix in VL
//      (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
//
        orghr(iLo, iHi, VL, tau, hrdWork);
//
//      Perform QR iteration, accumulating Schur vectors in VL
//      (Workspace: need N+1, prefer N+HSWORK (see comments) )
//
        iWork = iTau;
        auto hseqrWork = work(_(iWork, lWork));
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
//      (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
//
        orghr(iLo, iHi, VR, tau, hrdWork);
//
//      Perform QR iteration, accumulating Schur vectors in VR
//      (Workspace: need N+1, prefer N+HSWORK (see comments) )
//
        iWork = iTau;
        auto hseqrWork = work(_(iWork, lWork));
        info = hseqr(HSEQR::Schur, HSEQR::NoInit, iLo, iHi, A,
                     wr, wi, VR, hseqrWork);

    } else {
//
//      Compute eigenvalues only
//      (Workspace: need N+1, prefer N+HSWORK (see comments) )
//
        iWork = iTau;
        auto hseqrWork = work(_(iWork, lWork));
        info = hseqr(HSEQR::Eigenvalues, HSEQR::No, iLo, iHi, A,
                     wr, wi, VR, hseqrWork);
    }

//
//  If INFO > 0 from DHSEQR, then quit
//
    if (info==0) {

        if (computeVL || computeVR) {
//
//          Compute left and/or right eigenvectors
//          (Workspace: need 4*N)
//
            // TODO: shouldn't that be "Workspace: need 3*N" ?? instead of 4*N?
            IndexType nOut;
            auto trevcWork = work(_(iWork, lWork));

            // TODO: find a way that we don't have to pass an empty "select"
            DenseVector<Array<bool> >   select;
            trevc(computeVL, computeVR, TREVC::Backtransform, select,
                  A, VL, VR, n, nOut, trevcWork);
        }

        if (computeVL) {
//
//          Undo balancing of left eigenvectors
//          (Workspace: need N)
//
            bak(BALANCE::Both, Left, iLo, iHi, balWork, VL);
//
//          Normalize left eigenvectors and make largest component real
//
            auto _work = work(_(iWork, iWork+n-1));
            for (IndexType i=1; i<=n; ++i) {
                if (wi(i)==Zero) {
                    VL(_,i) *= One / blas::nrm2(VL(_,i));
                } else if (wi(i)>Zero) {
                    T scale = One / lapy2(blas::nrm2(VL(_,i  )),
                                          blas::nrm2(VL(_,i+1)));
                    VL(_,i)   *= scale;
                    VL(_,i+1) *= scale;
                    for (IndexType k=1; k<=n; ++k) {
                        _work(k) = pow(VL(k,i), 2) + pow(VL(k,i+1), 2);
                    }
                    IndexType k = blas::iamax(_work);
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
//          (Workspace: need N)
//
            bak(BALANCE::Both, Right, iLo, iHi, balWork, VR);
//
//          Normalize right eigenvectors and make largest component real
//
            auto _work = work(_(iWork, iWork+n-1));
            for (IndexType i=1; i<=n; ++i) {
                if (wi(i)==Zero) {
                    VR(_,i) *= One / blas::nrm2(VR(_,i));
                } else if (wi(i)>Zero) {
                    T scale = One / lapy2(blas::nrm2(VR(_,i  )),
                                          blas::nrm2(VR(_,i+1)));
                    VR(_,i)   *= scale;
                    VR(_,i+1) *= scale;
                    for (IndexType k=1; k<=n; ++k) {
                        _work(k) = pow(VR(k,i), 2) + pow(VR(k,i+1), 2);
                    }
                    IndexType k = blas::iamax(_work);
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
        lascl(LASCL::FullMatrix, 0, 0, cScale, ANorm, wr(_(info+1,n)));
        lascl(LASCL::FullMatrix, 0, 0, cScale, ANorm, wi(_(info+1,n)));

        if (info>0) {
            lascl(LASCL::FullMatrix, 0, 0, cScale, ANorm, wr(_(1,iLo-1)));
            lascl(LASCL::FullMatrix, 0, 0, cScale, ANorm, wi(_(1,iLo-1)));
        }
    }

    work(1) = maxWork;
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)ev_wsq [worksize query, real variant] ---------------------------------

template <typename MA>
typename RestrictTo<IsNotComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq_impl(bool computeVL, bool computeVR, const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.numRows();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        if (computeVL || computeVR) {
            minWork = 4*n;
        } else {
            minWork = 3*n;
        }
    }

//
//  Get optimal workspace from external LAPACK
//
    T           DUMMY, WORK;
    IndexType   LDVL = computeVL ? max(A.numRows(), IndexType(1)) : 1,
                LDVR = computeVR ? max(A.numRows(), IndexType(1)) : 1;
    IndexType   LWORK = -1;
    cxxlapack::geev<IndexType>(computeVL ? 'V' : 'N',
                               computeVR ? 'V' : 'N',
                               A.numRows(),
                               &DUMMY,
                               A.leadingDimension(),
                               &DUMMY,
                               &DUMMY,
                               &DUMMY,
                               LDVL,
                               &DUMMY,
                               LDVR,
                               &WORK,
                               LWORK);
    return Pair<IndexType>(minWork, WORK);
}

//-- (ge)ev_wsq [worksize query, complex variant] ------------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq_impl(bool computeVL, bool computeVR, const GeMatrix<MA> &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     RT;
    typedef typename GeMatrix<MA>::IndexType            IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.numRows();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        minWork = 2*n;
    }

//
//  Get optimal workspace from external LAPACK
//
    T           DUMMY, WORK;
    RT          RDUMMY;
    IndexType   LDVL = computeVL ? max(A.numRows(), IndexType(1)) : 1,
                LDVR = computeVR ? max(A.numRows(), IndexType(1)) : 1;
    IndexType   LWORK = -1;
    cxxlapack::geev<IndexType>(computeVL ? 'V' : 'N',
                               computeVR ? 'V' : 'N',
                               A.numRows(),
                               &DUMMY,
                               A.leadingDimension(),
                               &DUMMY,
                               &DUMMY,
                               LDVL,
                               &DUMMY,
                               LDVR,
                               &WORK,
                               LWORK,
                               &RDUMMY);
    return Pair<IndexType>(minWork, WORK.real());
}

//-- (ge)ev [real variant] -----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
typename GeMatrix<MA>::IndexType
ev_impl(bool                  computeVL,
        bool                  computeVR,
        GeMatrix<MA>          &A,
        DenseVector<VWR>      &wr,
        DenseVector<VWI>      &wi,
        GeMatrix<MVL>         &VL,
        GeMatrix<MVR>         &VR,
        DenseVector<VWORK>    &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    if (work.length()==0) {
        const auto ws = ev_wsq_impl(computeVL, computeVR, A);
        work.resize(ws.second, 1);
    }

    IndexType  info;
    info = cxxlapack::geev(computeVL ? 'V' : 'N',
                           computeVR ? 'V' : 'N',
                           A.numRows(),
                           A.data(),
                           A.leadingDimension(),
                           wr.data(),
                           wi.data(),
                           VL.data(),
                           VL.leadingDimension(),
                           VR.data(),
                           VR.leadingDimension(),
                           work.data(),
                           work.length());
    ASSERT(info>=0);
    return info;
}

//-- (ge)ev [complex variant] --------------------------------------------------

template <typename MA, typename VW, typename MVL, typename MVR, typename VWORK,
          typename VRWORK>
typename GeMatrix<MA>::IndexType
ev_impl(bool                  computeVL,
        bool                  computeVR,
        GeMatrix<MA>          &A,
        DenseVector<VW>       &w,
        GeMatrix<MVL>         &VL,
        GeMatrix<MVR>         &VR,
        DenseVector<VWORK>    &work,
        DenseVector<VRWORK>   &rWork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    if (work.length()==0) {
        const auto ws = ev_wsq_impl(computeVL, computeVR, A);
        work.resize(ws.second, 1);
    }

    IndexType  info;
    info = cxxlapack::geev(computeVL ? 'V' : 'N',
                           computeVR ? 'V' : 'N',
                           A.numRows(),
                           A.data(),
                           A.leadingDimension(),
                           w.data(),
                           VL.data(),
                           VL.leadingDimension(),
                           VR.data(),
                           VR.leadingDimension(),
                           work.data(),
                           work.length(),
                           rWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)ev [real variant] -----------------------------------------------------

template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MVL>::value
                 && IsRealGeMatrix<MVR>::value
                 && IsRealDenseVector<VWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeVL,
   bool     computeVR,
   MA       &&A,
   VWR      &&wr,
   VWI      &&wi,
   MVL      &&VL,
   MVR      &&VR,
   VWORK    &&work)
{
    LAPACK_DEBUG_OUT("(ge)ev [real]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<VWR>::Type    VectorWR;
    typedef typename RemoveRef<VWI>::Type    VectorWI;
    typedef typename RemoveRef<MVL>::Type    MatrixVL;
    typedef typename RemoveRef<MVR>::Type    MatrixVR;
    typedef typename RemoveRef<VWORK>::Type  VectorWork;

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

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView      A_org = A;
    typename VectorWR::NoView     wr_org   = wr;
    typename VectorWI::NoView     wi_org   = wi;
    typename MatrixVL::NoView     VL_org   = VL;
    typename MatrixVR::NoView     VR_org   = VR;
    typename VectorWork::NoView   work_org = work;
#   endif

//
//  Call implementation
//
    IndexType result = LAPACK_SELECT::ev_impl(computeVL, computeVR,
                                              A, wr, wi, VL, VR,
                                              work);
#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView      A_generic = A;
    typename VectorWR::NoView     wr_generic   = wr;
    typename VectorWI::NoView     wi_generic   = wi;
    typename MatrixVL::NoView     VL_generic   = VL;
    typename MatrixVR::NoView     VR_generic   = VR;
    typename VectorWork::NoView   work_generic = work;

    A    = A_org;
    wr   = wr_org;
    wi   = wi_org;
    VL   = VL_org;
    VR   = VR_org;
    work = work_org;

    IndexType _result = external::ev_impl(computeVL, computeVR,
                                          A, wr, wi, VL, VR,
                                          work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(wr_generic, wr, " wr_generic", "wr")) {
        std::cerr << "CXXLAPACK: wr_generic = " << wr_generic << std::endl;
        std::cerr << "F77LAPACK: wr = " << wr << std::endl;
        failed = true;
    }

    if (! isIdentical(wi_generic, wi, " wi_generic", "wi")) {
        std::cerr << "CXXLAPACK: wi_generic = " << wi_generic << std::endl;
        std::cerr << "F77LAPACK: wi = " << wi << std::endl;
        failed = true;
    }

    if (! isIdentical(VL_generic, VL, " VL_generic", "VL")) {
        std::cerr << "CXXLAPACK: VL_generic = " << VL_generic << std::endl;
        std::cerr << "F77LAPACK: VL = " << VL << std::endl;
        failed = true;
    }

    if (! isIdentical(VR_generic, VR, " VR_generic", "VR")) {
        std::cerr << "CXXLAPACK: VR_generic = " << VR_generic << std::endl;
        std::cerr << "F77LAPACK: VR = " << VR << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(result, _result, " result", "_result")) {
        std::cerr << "CXXLAPACK:  result = " << result << std::endl;
        std::cerr << "F77LAPACK: _result = " << _result << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
//        std::cerr << "passed: ev.tcc" << std::endl;
    }

#   endif

    return result;
}

//-- (ge)ev [complex variant] -----------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename VW, typename MVL, typename MVR, typename VWORK,
          typename VRWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MVL>::value
                 && IsComplexGeMatrix<MVR>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeVL,
   bool     computeVR,
   MA       &&A,
   VW       &&w,
   MVL      &&VL,
   MVR      &&VR,
   VWORK    &&work,
   VRWORK   &&rWork)
{
    LAPACK_DEBUG_OUT("(ge)ev [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;
    typedef typename RemoveRef<VW>::Type      VectorWR;
    typedef typename RemoveRef<MVL>::Type     MatrixVL;
    typedef typename RemoveRef<MVR>::Type     MatrixVR;
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
    typedef typename RemoveRef<VRWORK>::Type  VectorRWork;

    const IndexType n = A.numRows();

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rWork.firstIndex()==1);

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
#   endif

//
//  Resize output arguments if they are empty and needed
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

//
//  Call external implementation
//
    IndexType result = external::ev_impl(computeVL, computeVR, A, w, VL, VR,
                                         work, rWork);
    return result;
}

#endif // USE_CXXLAPACK

//-- (ge)ev_wsq [worksize query, real variant] ---------------------------------

template <typename MA>
typename RestrictTo<IsRealGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq(bool computeVL, bool computeVR, const MA &A)
{
    LAPACK_DEBUG_OUT("ev_wsq [real]");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif

//
//  Call implementation
//
    const auto ws = LAPACK_SELECT::ev_wsq_impl(computeVL, computeVR, A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    auto optWorkSize = external::ev_wsq_impl(computeVL, computeVR, A);
    if (! isIdentical(optWorkSize.second, ws.second,
                      "optWorkSize", "ws.second"))
    {
        ASSERT(0);
    }
#   endif

    return ws;
}

//-- (ge)ev_wsq [worksize query, complex variant] ------------------------------


#ifdef USE_CXXLAPACK

template <typename MA>
typename RestrictTo<IsComplexGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq(bool computeVL, bool computeVR, const MA &A)
{
    LAPACK_DEBUG_OUT("ev_wsq [complex]");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif

//
//  Call implementation
//
    const auto ws = external::ev_wsq_impl(computeVL, computeVR, A);

    return ws;
}

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_EV_TCC
