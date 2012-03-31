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
       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
      $                  VS, LDVS, WORK, LWORK, BWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_EIG_ES_TCC
#define FLENS_LAPACK_EIG_ES_TCC 1

#include <cmath>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA>
Pair<typename GeMatrix<MA>::IndexType>
es_generic_wsq(bool                 computeSchurVectors,
               const GeMatrix<MA>   &A)
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
        minWork = 3*n;

        HSEQR::ComputeZ computeZ = (computeSchurVectors) ? HSEQR::NoInit
                                                         : HSEQR::No;
        IndexType hsWork = hseqr_wsq(HSEQR::Schur, computeZ,
                                     IndexType(1), n, A);

        if (!computeSchurVectors) {
            maxWork = max(maxWork, n + hsWork);
        } else {
            maxWork = max(maxWork,
                          2*n + (n-1)*ilaenv<T>(1, "ORGHR", "", n, 1, n));
            maxWork = max(maxWork, n +hsWork);
        }
    }
    return Pair<IndexType>(minWork, maxWork);
}

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es_generic(bool                 computeSchurVectors,
           bool                 sortEigenvalues,
           SelectFunction       selectFunction,
           GeMatrix<MA>         &A,
           IndexType            &sDim,
           DenseVector<VWR>     &wr,
           DenseVector<VWI>     &wi,
           GeMatrix<MVS>        &VS,
           DenseVector<VWORK>   &work,
           DenseVector<BWORK>   &bWork)
{
    using std::sqrt;

    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;
    const IndexType n = A.numRows();
    const T Zero(0), One(1);

    IndexType info = 0;
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
    Pair<IndexType> wsQuery = es_wsq(computeSchurVectors, A);
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
        sDim = 0;
        return info;
    }
//
//  Get machine constants
//
    const T eps = lamch<T>(Precision);
    T smallNum  = lamch<T>(SafeMin);
    T bigNum    = One / smallNum;
    labad(smallNum, bigNum);

    smallNum    = sqrt(smallNum) / eps;
    bigNum      = One / smallNum;
//
//  Scale A if max element outside range [SMLNUM,BIGNUM]
//
    const T normA = lan(MaximumNorm, A);

    bool scaleA = false;
    T    cScale;

    if (normA>Zero && normA<smallNum) {
        scaleA = true;
        cScale = smallNum;
    } else if (normA>bigNum) {
        scaleA = true;
        cScale = bigNum;
    }
    if (scaleA) {
        lascl(LASCL::FullMatrix, 0, 0, normA, cScale, A);
    }
//
//  Permute the matrix to make it more nearly triangular
//  (Workspace: need N)
//
    IndexType iBal = 1;
    auto balWork = work(_(iBal,iBal+n-1));

    IndexType iLo, iHi;
    bal(BALANCE::PermuteOnly, A, iLo, iHi, balWork);
//
//  Reduce to upper Hessenberg form
//  (Workspace: need 3*N, prefer 2*N+N*NB)
//
    IndexType iTau  = iBal + n;
    IndexType iWork = iTau + n;
    IndexType lWork = work.length();

    auto tau = work(_(iTau, iTau+n-2));
    auto hrdWork = work(_(iWork, lWork));

    hrd(iLo, iHi, A, tau, hrdWork);

    if (computeSchurVectors) {
//
//      Copy Householder vectors to VS
//
        VS.lower() = A.lower();
//
//      Generate orthogonal matrix in VS
//      (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
//
        orghr(iLo, iHi, VS, tau, hrdWork);
    }

    sDim = 0;
//
//  Perform QR iteration, accumulating Schur vectors in VS if desired
//  (Workspace: need N+1, prefer N+HSWORK (see comments) )
//
    iWork = iTau;
    auto hseqrWork = work(_(iWork, lWork));
    HSEQR::ComputeZ computeZ = (computeSchurVectors) ? HSEQR::NoInit
                                                     : HSEQR::No;
    IndexType iEval = hseqr(HSEQR::Schur, computeZ, iLo, iHi,
                            A, wr, wi, VS, hseqrWork);
    if (iEval>0) {
        info = iEval;
    }
//
//  Sort eigenvalues if desired
//
    if (sortEigenvalues && info==0) {
        if (scaleA) {
            lascl(LASCL::FullMatrix, 0, 0, cScale, normA, wr);
            lascl(LASCL::FullMatrix, 0, 0, cScale, normA, wi);
        }
        for (IndexType i=1; i<=n; ++i) {
            bWork(i) = selectFunction(wr(i), wi(i));
        }
//
//      Reorder eigenvalues and transform Schur vectors
//      (Workspace: none needed)
//
        T sep, s;
        // TODO: I dislike that a dummy vector is needed
        DenseVector<Array<int> > dummy(1);
        IndexType iCond = trsen(TRSEN::None, computeSchurVectors, bWork,
                                A, VS, wr, wi, sDim, s, sep,
                                hseqrWork, dummy);
        if (iCond>0) {
            info = n + iCond;
        }
    }

    if (computeSchurVectors) {
//
//      Undo balancing
//      (Workspace: need N)
//
        bak(BALANCE::PermuteOnly, Right, iLo, iHi, balWork, VS);
    }

    if (scaleA) {
//
//      Undo scaling for the Schur form of A
//
        lascl(LASCL::UpperHessenberg, 0, 0, cScale, normA, A);
        wr = A.diag(0);
        if (cScale==smallNum) {
//
//          If scaling back towards underflow, adjust WI if an
//          offdiagonal element of a 2-by-2 block in the Schur form
//          underflows.
//
            IndexType i1, i2;
            if (iEval>0) {
                i1 = iEval + 1;
                i2 = iHi - 1;
                lascl(LASCL::FullMatrix, 0, 0, cScale, normA, wi(_(1,iLo-1)));
            } else if (sortEigenvalues) {
                i1 = 1;
                i2 = n - 1;
            } else {
                i1 = iLo;
                i2 = iHi - 1;
            }
            IndexType iNext = i1-1;
            for (IndexType i=i1; i<=i2; ++i) {
                if (i<iNext) {
                    continue;
                }
                if (wi(i)==Zero) {
                    iNext = i+1;
                } else {
                    if (A(i+1,i)==Zero) {
                        wi(i)   = Zero;
                        wi(i+1) = Zero;
                    } else if (A(i+1,i)!=Zero && A(i,i+1)==Zero) {
                        wi(i)   = Zero;
                        wi(i+1) = Zero;
                        if (i>1) {
                            blas::swap(A(_(1,i-1),i), A(_(1,i-1),i+1)); 
                        }
                        if (n>i+1) {
                            blas::swap(A(i,_(i+2,n)), A(i+1,_(i+2,n)));
                        }
                        if (computeSchurVectors) {
                            blas::swap(VS(_,i), VS(_,i+1));
                        }
                        A(i,i+1) = A(i+1,i);
                        A(i+1,i) = Zero;
                    }
                    iNext = i + 2;
                }
            }
        }
//
//      Undo scaling for the imaginary part of the eigenvalues
//
        lascl(LASCL::FullMatrix, 0, 0, cScale, normA, wi(_(iEval+1,n)));
    }

    if (sortEigenvalues && info==0) {
//
//      Check if reordering successful
//
        bool lastSelect = true;
        bool last2Select = true;
        sDim = 0;
        IndexType ip = 0;
        for (IndexType i=1; i<=n; ++i) {
            bool currentSelect = selectFunction(wr(i), wi(i));
            if (wi(i)==Zero) {
                if (currentSelect) {
                    ++sDim;
                }
                ip = 0;
                if (currentSelect && !lastSelect) {
                    info = n + 2;
                }
            } else {
                if (ip==1) {
//
//                  Last eigenvalue of conjugate pair
//
                    currentSelect = currentSelect || lastSelect;
                    lastSelect = currentSelect;
                    if (currentSelect) {
                        sDim += 2;
                    }
                    ip = -1;
                    if (currentSelect && !last2Select) {
                        info = n + 2;
                    }
                } else {
//
//                  First eigenvalue of conjugate pair
//
                    ip = 1;
                }
            }
            last2Select = lastSelect;
            lastSelect = currentSelect;
        }
    }

    work(1) = maxWork;
    return info;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA>
typename GeMatrix<MA>::IndexType
es_wsq(bool                 computeSchurVectors,
       const GeMatrix<MA>   &A)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  We need a bunch of dummy arguments ...
//
    int        (*SELECT)(const T*, const T*) = 0;
    IndexType   SDIM    = 0;
    T           DUMMY   = 0;
    IndexType   BWORK;
    IndexType   LDVS    = computeSchurVectors ? A.numRows() : 1;
//
//  ... and some arguments for the workspace query
//
    T                    WORK    = 0;
    const IndexType      LWORK   = -1;

    cxxlapack::gees<IndexType>(computeSchurVectors ? 'V' : 'N',
                               'N',
                               SELECT,
                               A.numRows(),
                               &DUMMY,
                               A.leadingDimension(),
                               SDIM,
                               &DUMMY,
                               &DUMMY,
                               &DUMMY,
                               LDVS,
                               &WORK,
                               LWORK,
                               &BWORK);
    return WORK;
}

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es(bool                 computeSchurVectors,
   bool                 sortEigenvalues,
   SelectFunction       selectFunction,
   GeMatrix<MA>         &A,
   IndexType            &sDim,
   DenseVector<VWR>     &wr,
   DenseVector<VWI>     &wi,
   GeMatrix<MVS>        &VS,
   DenseVector<VWORK>   &work,
   DenseVector<BWORK>   &bWork)
{
    IndexType info;
    info = cxxlapack::gees<IndexType>(computeSchurVectors ? 'V' : 'N',
                                      sortEigenvalues ? 'S' : 'N',
                                      selectFunction.select,
                                      A.numRows(),
                                      A.data(),
                                      A.leadingDimension(),
                                      sDim,
                                      wr.data(),
                                      wi.data(),
                                      VS.data(),
                                      VS.leadingDimension(),
                                      work.data(),
                                      work.length(),
                                      bWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA>
Pair<typename GeMatrix<MA>::IndexType>
es_wsq(bool                 computeSchurVectors,
       const GeMatrix<MA>   &A)
{
    LAPACK_DEBUG_OUT("es_wsq");

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
    const auto ws = es_generic_wsq(computeSchurVectors, A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    auto optWorkSize = external::es_wsq(computeSchurVectors, A);

    if (! isIdentical(optWorkSize, ws.second, "optWorkSize", "ws.second")) {
        ASSERT(0);
    }
#   endif

    return ws;

}

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es(bool                 computeSchurVectors,
   bool                 sortEigenvalues,
   SelectFunction       selectFunction,
   GeMatrix<MA>         &A,
   IndexType            &sDim,
   DenseVector<VWR>     &wr,
   DenseVector<VWI>     &wi,
   GeMatrix<MVS>        &VS,
   DenseVector<VWORK>   &work,
   DenseVector<BWORK>   &bWork)
{
    LAPACK_DEBUG_OUT("es");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);

    typename GeMatrix<MA>::IndexType n = A.numRows();

    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==n);

    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==n);

    if (computeSchurVectors) {
        ASSERT(VS.numRows()==VS.numCols());
        ASSERT(VS.numRows()==n);
        ASSERT(VS.firstRow()==1);
        ASSERT(VS.firstCol()==1);
    }

    if (sortEigenvalues) {
        ASSERT(bWork.firstIndex()==1);
        ASSERT(bWork.length()>=n);
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

    typename GeMatrix<MA>::NoView         A_org     = A;
    IndexType                             sDim_org  = sDim;
    typename DenseVector<VWR>::NoView     wr_org    = wr;
    typename DenseVector<VWI>::NoView     wi_org    = wi;
    typename GeMatrix<MVS>::NoView        VS_org    = VS;
    typename DenseVector<VWORK>::NoView   work_org  = work;
    typename DenseVector<BWORK>::NoView   bWork_org = bWork;

#   endif

//
//  Call implementation
//
    IndexType result = es_generic(computeSchurVectors,
                                  sortEigenvalues,
                                  selectFunction,
                                  A,
                                  sDim,
                                  wr,
                                  wi,
                                  VS,
                                  work,
                                  bWork);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView         A_generic     = A;
    IndexType                             sDim_generic  = sDim;
    typename DenseVector<VWR>::NoView     wr_generic    = wr;
    typename DenseVector<VWI>::NoView     wi_generic    = wi;
    typename GeMatrix<MVS>::NoView        VS_generic    = VS;
    typename DenseVector<VWORK>::NoView   work_generic  = work;
    typename DenseVector<BWORK>::NoView   bWork_generic = bWork;
//
//  restore output arguments
//
    A     = A_org;
    sDim  = sDim_org;
    wr    = wr_org;
    wi    = wi_org;
    VS    = VS_org;
    work  = work_org;
    bWork = bWork_org;

//
//  Compare generic results with results from the native implementation
//
    IndexType _result = external::es(computeSchurVectors,
                                     sortEigenvalues,
                                     selectFunction,
                                     A,
                                     sDim,
                                     wr,
                                     wi,
                                     VS,
                                     work,
                                     bWork);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(sDim_generic, sDim, "sDim_generic", "sDim")) {
        std::cerr << "CXXLAPACK: sDim_generic = " << sDim_generic << std::endl;
        std::cerr << "F77LAPACK: sDim = " << sDim << std::endl;
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

    if (! isIdentical(VS_generic, VS, "VS_generic", "VS")) {
        std::cerr << "CXXLAPACK: VS_generic = " << VS_generic << std::endl;
        std::cerr << "F77LAPACK: VS = " << VS << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(bWork_generic, bWork, "bWork_generic", "bWork")) {
        std::cerr << "CXXLAPACK: bWork_generic = "
                  << bWork_generic << std::endl;
        std::cerr << "F77LAPACK: bWork = " << bWork << std::endl;
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
//        std::cerr << "passed: es.tcc" << std::endl;
    }

#   endif

    return result;
}

//-- forwarding ----------------------------------------------------------------
template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es(bool                 computeSchurVectors,
   bool                 sortEigenvalues,
   SelectFunction       selectFunction,
   MA                   &&A,
   IndexType            &&sDim,
   VWR                  &&wr,
   VWI                  &&wi,
   MVS                  &&VS,
   VWORK                &&work,
   BWORK                &&bWork)
{
    CHECKPOINT_ENTER;
    const IndexType info = es(computeSchurVectors,
                              sortEigenvalues,
                              selectFunction,
                              A,
                              sDim,
                              wr,
                              wi,
                              VS,
                              work,
                              bWork);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_ES_TCC
