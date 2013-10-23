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
      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
     $                   LDZ, WORK, LWORK, INFO )
      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
     $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK computational routine (version 3.2.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     June 2010
 */

#ifndef FLENS_LAPACK_IMPL_HSEQR_TCC
#define FLENS_LAPACK_IMPL_HSEQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename T>
void
_makeHseqrOpt(char opt[3], HSEQR::Job job, HSEQR::ComputeZ computeZ)
{
    opt[0] = (job==HSEQR::Eigenvalues) ? 'E' : 'S';
    if (computeZ==HSEQR::No) {
        opt[1] = 'N';
    } else if (computeZ==HSEQR::Init) {
        opt[1] = 'I';
    } else {
        opt[1] = 'V';
    }
    opt[2] = 0;
}

namespace generic {

//
//  Workspace query (real/complex variant)
//
template <typename IndexType, typename MH>
IndexType
hseqr_wsq_impl(HSEQR::Job            job,
               HSEQR::ComputeZ       computeZ,
               IndexType             iLo,
               IndexType             iHi,
               const GeMatrix<MH>    &H)
{
    using std::max;

    IndexType n = H.numRows();
    const bool wantT = (job==HSEQR::Schur);
    const bool wantZ = (computeZ!=HSEQR::No);

    IndexType info = laqr0_wsq(wantT, wantZ, iLo, iHi, H);
    info = max(max(IndexType(1), n), info);
    return info;
}

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
hseqr_impl(HSEQR::Job            job,
           HSEQR::ComputeZ       computeZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VWR>      &wr,
           DenseVector<VWI>      &wi,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType   T;
    typedef typename GeMatrix<MH>::View          GeMatrixView;
    typedef typename GeMatrix<MH>::VectorView    DenseVectorView;

    const IndexType nTiny = 11;
//  ==== NL allocates some local workspace to help small matrices
//  .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
//  .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
//  .    mended.  (The default value of NMIN is 75.)  Using NL = 49
//  .    allows up to six simultaneous shifts and a 16-by-16
//  .    deflation window.  ====
    const           IndexType nl = 49;
    T               hlBuffer[nl*nl], worklBuffer[nl];
    GeMatrixView    Hl    = typename GeMatrixView::Engine(nl, nl, hlBuffer, nl);
    DenseVectorView workl = typename DenseVectorView::Engine(nl, worklBuffer);

    const Underscore<IndexType>     _;
    const IndexType                 n = H.numCols();
    const T                         Zero(0), One(1);

    IndexType info = 0;

//
//  ==== Decode and check the input parameters. ====
//
    const bool wantT = (job==HSEQR::Schur);
    const bool initZ = (computeZ==HSEQR::Init);
    const bool wantZ = (computeZ!=HSEQR::No);

    if (work.length()>0) {
        work(1) = n;
    } else {
//
//      ==== Perform and apply a workspace query ====
//
        IndexType lWorkOpt = laqr0_wsq(wantT, wantZ, iLo, iHi, H);
        work.resize(lWorkOpt);
    }

    if (n==0) {
//
//      ==== Quick return in case n = 0; nothing to do. ====
//
        return info;
    }
//
//  ==== copy eigenvalues isolated by DGEBAL ====
//
    for (IndexType i=1; i<=iLo-1; ++i) {
        wr(i) = H(i,i);
        wi(i) = Zero;
    }

    for (IndexType i=iHi+1; i<=n; ++i) {
        wr(i) = H(i,i);
        wi(i) = Zero;
    }
//
//  ==== Initialize Z, if requested ====
//
    if (initZ) {
        Z = Zero;
        Z.diag(0) = One;
    }
//
//  ==== Quick return if possible ====
//
    if (iLo==iHi) {
        wr(iLo) = H(iLo, iLo);
        wi(iLo) = Zero;
        return info;
    }
//
//  ==== lahqr/laqr0 crossover point ====
//
    char opt[3];
    _makeHseqrOpt<T>(opt, job, computeZ);
    IndexType nMin = ilaenv<T>(12, "HSEQR", opt, n, iLo, iHi, work.length());

    nMin= max(nTiny, nMin);
//
//  ==== laqr0 for big matrices; lahqr for small ones ====
//
    if (n>nMin) {
        info = laqr0(wantT, wantZ, iLo, iHi, H, wr, wi, iLo, iHi, Z, work);
    } else {
//
//      ==== Small matrix ====
//
        info = lahqr(wantT, wantZ, iLo, iHi, H, wr, wi, iLo, iHi, Z);

        if (info>0) {
//
//          ==== A rare lahqr failure!  laqr0 sometimes succeeds
//          .    when lahqr fails. ====
//
            const IndexType kBot = info;
            if (n>=nl) {
//
//              ==== Larger matrices have enough subdiagonal scratch
//              .    space to call laqr0 directly. ====
//
                info = laqr0(wantT, wantZ,
                             iLo, kBot, H, wr, wi,
                             iLo, iHi, Z, work);
            } else {
//
//              ==== Tiny matrices don't have enough subdiagonal
//              .    scratch space to benefit from DLAQR0.  Hence,
//              .    tiny matrices must be copied into a larger
//              .    array before calling DLAQR0. ====
//
                Hl(_(1,n),_(1,n)) = H;
                Hl(n+1, n) = Zero;
                Hl(_(1,nl),_(n+1,nl)) = Zero;
                info = laqr0(wantT, wantZ,
                             iLo, kBot, Hl, wr, wi,
                             iLo, iHi, Z, workl);
                if (wantT || (info!=0)) {
                    H = Hl(_(1,n),_(1,n));
                }
            }
        }
    }
//
//  ==== Clear out the trash, if necessary. ====
//
    if ((wantT || (info!=0)) && (n>2)) {
        H(_(3,n),_(1,n-2)).lower() = Zero;
    }
//
//  ==== Ensure reported workspace size is backward-compatible with
//  .    previous LAPACK versions. ====
//
    work(1) = max(T(max(IndexType(1),n)), work(1));
    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
IndexType
hseqr_impl(HSEQR::Job            job,
           HSEQR::ComputeZ       computeZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VW>       &w,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef typename GeMatrix<MH>::View                 GeMatrixView;
    typedef typename GeMatrix<MH>::VectorView           DenseVectorView;

    const IndexType nTiny = 11;
//  ==== NL allocates some local workspace to help small matrices
//  .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
//  .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
//  .    mended.  (The default value of NMIN is 75.)  Using NL = 49
//  .    allows up to six simultaneous shifts and a 16-by-16
//  .    deflation window.  ====
    const           IndexType nl = 49;
    T               hlBuffer[nl*nl], worklBuffer[nl];
    GeMatrixView    Hl    = typename GeMatrixView::Engine(nl, nl, hlBuffer, nl);
    DenseVectorView workl = typename DenseVectorView::Engine(nl, worklBuffer);

    const Underscore<IndexType>     _;
    const IndexType                 n = H.numCols();
    const T                         Zero(0), One(1);

    IndexType info = 0;

//
//  ==== Decode and check the input parameters. ====
//
    const bool wantT = (job==HSEQR::Schur);
    const bool initZ = (computeZ==HSEQR::Init);
    const bool wantZ = (computeZ!=HSEQR::No);

    if (work.length()>0) {
        work(1) = n;
    } else {
//
//      ==== Perform and apply a workspace query ====
//
        IndexType lWorkOpt = laqr0_wsq(wantT, wantZ, iLo, iHi, H);
        work.resize(lWorkOpt);
    }

    if (n==0) {
//
//      ==== Quick return in case n = 0; nothing to do. ====
//
        return info;
    }
//
//  ==== copy eigenvalues isolated by DGEBAL ====
//
    for (IndexType i=1; i<=iLo-1; ++i) {
        w(i) = H(i,i);
    }

    for (IndexType i=iHi+1; i<=n; ++i) {
        w(i) = H(i,i);
    }
//
//  ==== Initialize Z, if requested ====
//
    if (initZ) {
        Z = Zero;
        Z.diag(0) = One;
    }
//
//  ==== Quick return if possible ====
//
    if (iLo==iHi) {
        w(iLo) = H(iLo, iLo);
        return info;
    }
//
//  ==== lahqr/laqr0 crossover point ====
//
    char opt[3];
    _makeHseqrOpt<T>(opt, job, computeZ);
    IndexType nMin = ilaenv<T>(12, "HSEQR", opt, n, iLo, iHi, work.length());

    nMin= max(nTiny, nMin);
//
//  ==== laqr0 for big matrices; lahqr for small ones ====
//
    if (n>nMin) {
        info = laqr0(wantT, wantZ, iLo, iHi, H, w, iLo, iHi, Z, work);
    } else {
//
//      ==== Small matrix ====
//
        info = lahqr(wantT, wantZ, iLo, iHi, H, w, iLo, iHi, Z);

        if (info>0) {
//
//          ==== A rare lahqr failure!  laqr0 sometimes succeeds
//          .    when lahqr fails. ====
//
            const IndexType kBot = info;
            if (n>=nl) {
//
//              ==== Larger matrices have enough subdiagonal scratch
//              .    space to call laqr0 directly. ====
//
                info = laqr0(wantT, wantZ,
                             iLo, kBot, H, w,
                             iLo, iHi, Z, work);
            } else {
//
//              ==== Tiny matrices don't have enough subdiagonal
//              .    scratch space to benefit from DLAQR0.  Hence,
//              .    tiny matrices must be copied into a larger
//              .    array before calling DLAQR0. ====
//
                Hl(_(1,n),_(1,n)) = H;
                Hl(n+1, n) = Zero;
                Hl(_(1,nl),_(n+1,nl)) = Zero;
                info = laqr0(wantT, wantZ,
                             iLo, kBot, Hl, w,
                             iLo, iHi, Z, workl);
                if (wantT || (info!=0)) {
                    H = Hl(_(1,n),_(1,n));
                }
            }
        }
    }
//
//  ==== Clear out the trash, if necessary. ====
//
    if ((wantT || (info!=0)) && (n>2)) {
        H(_(3,n),_(1,n-2)).lower() = Zero;
    }
//
//  ==== Ensure reported workspace size is backward-compatible with
//  .    previous LAPACK versions. ====
//
    work(1) = T(max(PT(max(IndexType(1),n)),
                    cxxblas::real(work(1))));
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Workspace query (real variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsReal<typename MH::ElementType>::value,
         IndexType>::Type
hseqr_wsq_impl(HSEQR::Job         job,
               HSEQR::ComputeZ    computeZ,
               IndexType          iLo,
               IndexType          iHi,
               const GeMatrix<MH> &H)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    T               WORK  = 0;
    T               DUMMY = 0;
    const IndexType LWORK = -1;

    cxxlapack::hseqr<IndexType>(getF77Char(job),
                                getF77Char(computeZ),
                                H.numCols(),
                                iLo,
                                iHi,
                                &DUMMY,
                                H.leadingDimension(),
                                &DUMMY,
                                &DUMMY,
                                &DUMMY,
                                H.leadingDimension(),
                                &WORK,
                                LWORK);
    return WORK;
}

//
//  Workspace query (complex variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsComplex<typename MH::ElementType>::value,
         IndexType>::Type
hseqr_wsq_impl(HSEQR::Job         job,
               HSEQR::ComputeZ    computeZ,
               IndexType          iLo,
               IndexType          iHi,
               const GeMatrix<MH> &H)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    T               WORK  = 0;
    T               DUMMY = 0;
    const IndexType LWORK = -1;

    cxxlapack::hseqr<IndexType>(getF77Char(job),
                                getF77Char(computeZ),
                                H.numCols(),
                                iLo,
                                iHi,
                                &DUMMY,
                                H.leadingDimension(),
                                &DUMMY,
                                &DUMMY,
                                H.leadingDimension(),
                                &WORK,
                                LWORK);
    return cxxblas::real(WORK);
}

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
hseqr_impl(HSEQR::Job            job,
           HSEQR::ComputeZ       computeZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VWR>      &wr,
           DenseVector<VWI>      &wi,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    IndexType  info;
    info = cxxlapack::hseqr<IndexType>(getF77Char(job),
                                       getF77Char(computeZ),
                                       H.numCols(),
                                       iLo,
                                       iHi,
                                       H.data(),
                                       H.leadingDimension(),
                                       wr.data(),
                                       wi.data(),
                                       Z.data(),
                                       Z.leadingDimension(),
                                       work.data(),
                                       work.length());
    ASSERT(info>=0);
    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
IndexType
hseqr_impl(HSEQR::Job            job,
           HSEQR::ComputeZ       computeZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VW>       &w,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    IndexType  info;
    info = cxxlapack::hseqr<IndexType>(getF77Char(job),
                                       getF77Char(computeZ),
                                       H.numCols(),
                                       iLo,
                                       iHi,
                                       H.data(),
                                       H.leadingDimension(),
                                       w.data(),
                                       Z.data(),
                                       Z.leadingDimension(),
                                       work.data(),
                                       work.length());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//
//  Workspace query (real variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsRealGeMatrix<MH>::value,
         IndexType>::Type
hseqr_wsq(HSEQR::Job            job,
          HSEQR::ComputeZ       computeZ,
          IndexType             iLo,
          IndexType             iHi,
          const MH              &H)
{
    LAPACK_DEBUG_OUT("hseqr_wsq");

//
//  Test the input parameters
//
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::hseqr_wsq_impl(job, computeZ, iLo, iHi, H);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info =  external::hseqr_wsq_impl(job, computeZ, iLo, iHi, H);

    if (info!=_info) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        ASSERT(0);
    }
#   endif
    return info;
}

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MH>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MZ>::value
                 && IsRealDenseVector<VWORK>::value,
         IndexType>::Type
hseqr(HSEQR::Job                job,
      HSEQR::ComputeZ           computeZ,
      IndexType                 iLo,
      IndexType                 iHi,
      MH                        &&H,
      VWR                       &&wr,
      VWI                       &&wi,
      MZ                        &&Z,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("hseqr");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type     MatrixH;
    typedef typename RemoveRef<VWR>::Type    VectorWR;
    typedef typename RemoveRef<VWI>::Type    VectorWI;
    typedef typename RemoveRef<MZ>::Type     MatrixZ;
    typedef typename RemoveRef<VWORK>::Type  VectorWork;
#   endif

//
//  Test the input parameters
//
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());
    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==H.numCols());
    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==H.numCols());
    ASSERT((computeZ==HSEQR::No) || (Z.numRows()==H.numCols()));
    ASSERT((computeZ==HSEQR::No) || (Z.numCols()==H.numCols()));
    ASSERT((work.length()==0) || (work.length()>=H.numCols()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixH::NoView        H_org   = H;

    typename MatrixH::NoView        _H      = H;
    typename VectorWR::NoView       _wr     = wr;
    typename VectorWI::NoView       _wi     = wi;
    typename MatrixZ::NoView        _Z      = Z;
    typename VectorWork::NoView     _work   = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::hseqr_impl(job, computeZ, iLo, iHi, H,
                                               wr, wi, Z, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    // TODO: also check workspace query directly!
    if (_work.length()==0) {
        _work.resize(work.length());
    }

    IndexType _info = external::hseqr_impl(job, computeZ, iLo, iHi, _H,
                                           _wr, _wi, _Z, _work);

    bool failed = false;
    if (! isIdentical(H, _H, " H", "_H")) {
        std::cerr << "CXXLAPACK:  H = " << H << std::endl;
        std::cerr << "F77LAPACK: _H = " << _H << std::endl;
        failed = true;
    }

    if (! isIdentical(wr, _wr, " wr", "_wr")) {
        std::cerr << "CXXLAPACK:  wr = " << wr << std::endl;
        std::cerr << "F77LAPACK: _wr = " << _wr << std::endl;
        failed = true;
    }

    if (! isIdentical(wi, _wi, " wi", "_wi")) {
        std::cerr << "CXXLAPACK:  wi = " << wi << std::endl;
        std::cerr << "F77LAPACK: _wi = " << _wi << std::endl;
        failed = true;
    }

    if (! isIdentical(Z, _Z, " Z", "_Z")) {
        std::cerr << "CXXLAPACK:  Z = " << Z << std::endl;
        std::cerr << "F77LAPACK: _Z = " << _Z << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "H_org = " << H_org << std::endl;

        std::cerr << "job = " << job
                  << ", computeZ = " << computeZ
                  << ", iLo = " << iLo
                  << ", iHi = " << iHi
                  << ", H.numRows() = " << H.numRows()
                  << std::endl;
        ASSERT(0);
    }
#   endif
    return info;
}

//
//  Workspace query (complex variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsComplexGeMatrix<MH>::value,
         IndexType>::Type
hseqr_wsq(HSEQR::Job            job,
          HSEQR::ComputeZ       computeZ,
          IndexType             iLo,
          IndexType             iHi,
          const MH              &H)
{
    LAPACK_DEBUG_OUT("hseqr_wsq (complex)");

//
//  Test the input parameters
//
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::hseqr_wsq_impl(job, computeZ, iLo, iHi, H);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info =  external::hseqr_wsq_impl(job, computeZ, iLo, iHi, H);

    if (info!=_info) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        ASSERT(0);
    }
#   endif
    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MH>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MZ>::value
                 && IsComplexDenseVector<VWORK>::value,
         IndexType>::Type
hseqr(HSEQR::Job                job,
      HSEQR::ComputeZ           computeZ,
      IndexType                 iLo,
      IndexType                 iHi,
      MH                        &&H,
      VW                        &&w,
      MZ                        &&Z,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("hseqr (complex)");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type     MatrixH;
    typedef typename RemoveRef<VW>::Type     VectorW;
    typedef typename RemoveRef<MZ>::Type     MatrixZ;
    typedef typename RemoveRef<VWORK>::Type  VectorWork;
#   endif

//
//  Test the input parameters
//
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());
    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==H.numCols());
    ASSERT((computeZ==HSEQR::No) || (Z.numRows()==H.numCols()));
    ASSERT((computeZ==HSEQR::No) || (Z.numCols()==H.numCols()));
    ASSERT((work.length()==0) || (work.length()>=H.numCols()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixH::NoView        H_org   = H;

    typename MatrixH::NoView        _H      = H;
    typename VectorW::NoView        _w      = w;
    typename MatrixZ::NoView        _Z      = Z;
    typename VectorWork::NoView     _work   = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::hseqr_impl(job, computeZ, iLo, iHi, H,
                                               w, Z, work);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    // TODO: also check workspace query directly!
    if (_work.length()==0) {
        _work.resize(work.length());
    }

    IndexType _info = external::hseqr_impl(job, computeZ, iLo, iHi, _H,
                                           _w, _Z, _work);

    bool failed = false;
    if (! isIdentical(H, _H, " H", "_H")) {
        std::cerr << "CXXLAPACK:  H = " << H << std::endl;
        std::cerr << "F77LAPACK: _H = " << _H << std::endl;
        failed = true;
    }

    if (! isIdentical(w, _w, " w", "_w")) {
        std::cerr << "CXXLAPACK:  w = " << w << std::endl;
        std::cerr << "F77LAPACK: _w = " << _w << std::endl;
        failed = true;
    }

    if (! isIdentical(Z, _Z, " Z", "_Z")) {
        std::cerr << "CXXLAPACK:  Z = " << Z << std::endl;
        std::cerr << "F77LAPACK: _Z = " << _Z << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (! isIdentical(work, _work, " work", "_work")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: _work = " << _work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "H_org = " << H_org << std::endl;

        std::cerr << "job = " << job
                  << ", computeZ = " << computeZ
                  << ", iLo = " << iLo
                  << ", iHi = " << iHi
                  << ", H.numRows() = " << H.numRows()
                  << std::endl;
        ASSERT(0);
    }
#   endif
    return info;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_HSEQR_TCC
