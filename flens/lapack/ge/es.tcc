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
       SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
      $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_GE_ES_TCC
#define FLENS_LAPACK_GE_ES_TCC 1

#include <cxxstd/cmath.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (ge)es_wsq [real variant] -------------------------------------------------

template <typename MA>
Pair<typename GeMatrix<MA>::IndexType>
es_wsq_impl(bool                 computeSchurVectors,
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

//-- (ge)es [real variant] -----------------------------------------------------

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es_impl(bool                 computeSchurVectors,
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
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0), normA, cScale, A);
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
            lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                  cScale, normA, wr);
            lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                  cScale, normA, wi);
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
        DenseVector<Array<IndexType> >   dummy(1);
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
        lascl(LASCL::UpperHessenberg, IndexType(0), IndexType(0),
              cScale, normA, A);
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
                lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
                      cScale, normA, wi(_(1,iLo-1)));
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
        lascl(LASCL::FullMatrix, IndexType(0), IndexType(0),
              cScale, normA, wi(_(iEval+1,n)));
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

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (ge)es_wsq [real variant] -------------------------------------------------

template <typename MA>
typename RestrictTo<IsNotComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
es_wsq_impl(bool                 computeSchurVectors,
            const GeMatrix<MA>   &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.numRows();
    IndexType  minWork = (n==0) ? 1
                                : 3*n;

//
//  We need a bunch of dummy arguments ...
//
    IndexType   SDIM    = 0;
    T           DUMMY   = 0;
    IndexType   BWORK;
    IndexType   LDVS    = computeSchurVectors ? max(A.numRows(), IndexType(1))
                                              : 1;
//
//  ... and some arguments for the workspace query
//
    T                    WORK    = 0;
    const IndexType      LWORK   = -1;

    cxxlapack::gees<IndexType>(computeSchurVectors ? 'V' : 'N',
                               'N',
                               0,
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
    return Pair<IndexType>(minWork, WORK);
}

//-- (ge)es_wsq [complex variant] ----------------------------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
es_wsq_impl(bool                 computeSchurVectors,
            const GeMatrix<MA>   &A)
{
    using std::max;

    typedef typename GeMatrix<MA>::ElementType        T;
    typedef typename ComplexTrait<T>::PrimitiveType   RT;
    typedef typename GeMatrix<MA>::IndexType          IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.numRows();
    IndexType  minWork = (n==0) ? 1
                                : 3*n;

//
//  We need a bunch of dummy arguments ...
//
    IndexType   SDIM   = 0;
    T           DUMMY  = 0;
    RT          RDUMMY = 0;
    IndexType   BWORK;
    IndexType   LDVS    = computeSchurVectors ? max(A.numRows(), IndexType(1))
                                              : 1;
//
//  ... and some arguments for the workspace query
//
    T                    WORK    = 0;
    const IndexType      LWORK   = -1;

    cxxlapack::gees<IndexType>(computeSchurVectors ? 'V' : 'N',
                               'N',
                               0,
                               A.numRows(),
                               &DUMMY,
                               A.leadingDimension(),
                               SDIM,
                               &DUMMY,
                               &DUMMY,
                               LDVS,
                               &WORK,
                               LWORK,
                               &RDUMMY,
                               &BWORK);
    return Pair<IndexType>(minWork, WORK.real());
}

//-- (ge)es [real variant] -----------------------------------------------------

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename BWORK>
IndexType
es_impl(bool                 computeSchurVectors,
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
    if (work.length()==0) {
        const auto ws = es_wsq_impl(computeSchurVectors, A);
        work.resize(ws.second, 1);
    }

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

//-- (ge)es [complex variant] --------------------------------------------------

template <typename SelectFunction, typename MA, typename IndexType,
          typename VW, typename MVS, typename VWORK, typename VRWORK,
          typename VBWORK>
IndexType
es_impl(bool                 computeSchurVectors,
        bool                 sortEigenvalues,
        SelectFunction       selectFunction,
        GeMatrix<MA>         &A,
        IndexType            &sDim,
        DenseVector<VW>      &w,
        GeMatrix<MVS>        &VS,
        DenseVector<VWORK>   &work,
        DenseVector<VRWORK>  &rWork,
        DenseVector<VBWORK>  &bWork)
{
    if (work.length()==0) {
        const auto ws = es_wsq_impl(computeSchurVectors, A);
        work.resize(ws.second, 1);
    }

    IndexType info;
    info = cxxlapack::gees<IndexType>(computeSchurVectors ? 'V' : 'N',
                                      sortEigenvalues ? 'S' : 'N',
                                      selectFunction.select,
                                      A.numRows(),
                                      A.data(),
                                      A.leadingDimension(),
                                      sDim,
                                      w.data(),
                                      VS.data(),
                                      VS.leadingDimension(),
                                      work.data(),
                                      work.length(),
                                      rWork.data(),
                                      bWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (ge)es [real variant] -----------------------------------------------------

template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename VBWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsInteger<IndexType>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MVS>::value
                 && IsRealDenseVector<VWORK>::value
                 && IsRealDenseVector<VBWORK>::value,
         IndexType>::Type
es(bool                 computeSchurVectors,
   bool                 sortEigenvalues,
   SelectFunction       selectFunction,
   MA                   &&A,
   IndexType            &sDim,
   VWR                  &&wr,
   VWI                  &&wi,
   MVS                  &&VS,
   VWORK                &&work,
   VBWORK               &&bWork)
{
    LAPACK_DEBUG_OUT("(ge)es [real]");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type        MatrixA;
    typedef typename RemoveRef<VWR>::Type       VectorWR;
    typedef typename RemoveRef<VWI>::Type       VectorWI;
    typedef typename RemoveRef<MVS>::Type       MatrixVS;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
    typedef typename RemoveRef<VBWORK>::Type    VectorBWork;
#   endif

//
//  Test the input parameters
//
    const IndexType n = A.numRows();

#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(bWork.firstIndex()==1);

    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==0 || wr.length()==n);

    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==0 || wi.length()==n);

    if (computeSchurVectors) {
        ASSERT(VS.numRows()==VS.numCols());
        ASSERT(VS.numRows()==0 || VS.numRows()==n);
        ASSERT(VS.firstRow()==1);
        ASSERT(VS.firstCol()==1);
    }

    const auto ws = es_wsq(computeSchurVectors, A);
    ASSERT(work.length()==0 || work.length()>=ws.first);

    if (sortEigenvalues) {
        ASSERT(bWork.firstIndex()==1);
        ASSERT(bWork.length()==0 || bWork.length()>=n);
    }
#   endif

//
//  Resize output arguments.
//  Note: work gets resized by the actual implementation.
//
    if (wr.length()==0) {
        wr.resize(n, 1);
    }
    if (wi.length()==0) {
        wi.resize(n, 1);
    }
    if (computeSchurVectors && VS.numRows()==0) {
        VS.resize(n, n, 1, 1);
    }
    if (sortEigenvalues && bWork.length()==0) {
        bWork.resize(n, 1);
    }

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

    typename MatrixA::NoView      A_org     = A;
    IndexType                     sDim_org  = sDim;
    typename VectorWR::NoView     wr_org    = wr;
    typename VectorWI::NoView     wi_org    = wi;
    typename MatrixVS::NoView     VS_org    = VS;
    typename VectorWork::NoView   work_org  = work;
    typename VectorBWork::NoView  bWork_org = bWork;

#   endif

//
//  Call implementation
//
    IndexType result = LAPACK_SELECT::es_impl(computeSchurVectors,
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
    typename MatrixA::NoView      A_generic     = A;
    IndexType                     sDim_generic  = sDim;
    typename VectorWR::NoView     wr_generic    = wr;
    typename VectorWI::NoView     wi_generic    = wi;
    typename MatrixVS::NoView     VS_generic    = VS;
    typename VectorWork::NoView   work_generic  = work;
    typename VectorBWork::NoView  bWork_generic = bWork;
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
    IndexType result_ = external::es_impl(computeSchurVectors,
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

    if (! isIdentical(result, result_, " result", "result_")) {
        std::cerr << "CXXLAPACK:  result = " << result << std::endl;
        std::cerr << "F77LAPACK: result_ = " << result_ << std::endl;
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

#ifdef USE_CXXLAPACK

//-- (ge)es [complex variant] --------------------------------------------------

template <typename SelectFunction, typename MA, typename IndexType,
          typename VW, typename MVS, typename VWORK, typename VRWORK,
          typename VBWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsInteger<IndexType>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MVS>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value
                 && IsRealDenseVector<VBWORK>::value,
         IndexType>::Type
es(bool                 computeSchurVectors,
   bool                 sortEigenvalues,
   SelectFunction       select,
   MA                   &&A,
   IndexType            &sDim,
   VW                   &&w,
   MVS                  &&VS,
   VWORK                &&work,
   VRWORK               &&rWork,
   VBWORK               &&bWork)
{
    LAPACK_DEBUG_OUT("(ge)es [complex]");

//
//  Test the input parameters
//
    const IndexType n = A.numRows();

#   ifndef NDEBUG
    ASSERT(A.numRows()==A.numCols());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rWork.firstIndex()==1);
    ASSERT(bWork.firstIndex()==1);

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==0 || w.length()==n);

    if (computeSchurVectors) {
        ASSERT(VS.numRows()==VS.numCols());
        ASSERT(VS.numRows()==0 || VS.numRows()==n);
        ASSERT(VS.firstRow()==1);
        ASSERT(VS.firstCol()==1);
    }

    const auto ws = es_wsq(computeSchurVectors, A);
    ASSERT(work.length()==0 || work.length()>=ws.first);

    ASSERT(rWork.firstIndex()==1);
    ASSERT(rWork.length()==0 || rWork.length()==n);

    if (sortEigenvalues) {
        ASSERT(bWork.firstIndex()==1);
        ASSERT(bWork.length()==0 || bWork.length()>=n);
    }
#   endif

//
//  Resize output arguments.
//  Note: work gets resized by the actual implementation.
//
    if (w.length()==0) {
        w.resize(n, 1);
    }
    if (computeSchurVectors && VS.numRows()==0) {
        VS.resize(n, n, 1, 1);
    }
    if (rWork.length()==0) {
        rWork.resize(n, 1);
    }
    if (sortEigenvalues && bWork.length()==0) {
        bWork.resize(n, 1);
    }

    IndexType result_ = external::es_impl(computeSchurVectors,
                                          sortEigenvalues,
                                          select,
                                          A,
                                          sDim,
                                          w,
                                          VS,
                                          work,
                                          rWork,
                                          bWork);
    return result_;
}

#endif // USE_CXXLAPACK

//-- (ge)es_wsq [real variant] -------------------------------------------------

template <typename MA>
typename RestrictTo<IsRealGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
es_wsq(bool     computeSchurVectors,
       const MA &A)
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
    const auto ws = LAPACK_SELECT::es_wsq_impl(computeSchurVectors, A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const auto ws_ = external::es_wsq_impl(computeSchurVectors, A);

    if (! isIdentical(ws_.first, ws_.first, "ws_.first", "ws_.first")) {
        ASSERT(0);
    }
    if (! isIdentical(ws_.second, ws_.second, "ws_.second", "ws_.second")) {
        ASSERT(0);
    }
#   endif

    return ws;
}

//-- (ge)es_wsq [complex variant] ----------------------------------------------


#ifdef USE_CXXLAPACK

template <typename MA>
typename RestrictTo<IsComplexGeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
es_wsq(bool     computeSchurVectors,
       const MA &A)
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
    const auto ws = external::es_wsq_impl(computeSchurVectors, A);

    return ws;
}

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_ES_TCC
