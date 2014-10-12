/*
 *   Copyright (c) 2012, Klaus Pototzky
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
      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HE_EV_TCC
#define FLENS_LAPACK_HE_EV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (he)ev [worksize query hermitian variant] ---------------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq_impl(const HeMatrix<MA> &A)
{
    using std::max;

    typedef typename HeMatrix<MA>::ElementType  T;
    typedef typename HeMatrix<MA>::IndexType    IndexType;

    const IndexType n = A.numRows();

    const char      upLo[1]  = { getF77BlasChar(A.upLo()) };
    const IndexType nb       = ilaenv<T>(1, "HETRD", upLo, n);

    const IndexType lWorkOpt = max(IndexType(1), (nb+1)*n);
    const IndexType minWork  = (n==0) ? 1 : max(1,2*n-1);

    return Pair<IndexType>(minWork, lWorkOpt);
}

template <typename MA, typename VW, typename VWORK, typename VRWORK>
typename HeMatrix<MA>::IndexType
ev_impl(bool                  computeV,
        HeMatrix<MA>          &A,
        DenseVector<VW>       &w,
        DenseVector<VWORK>    &work,
        DenseVector<VRWORK>   &rWork)
{
    using std::real;

    typedef typename HeMatrix<MA>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef typename GeMatrix<MA>::IndexType            IndexType;

    const Underscore<IndexType> _;

    const bool upper  = (A.upLo()==Upper);

    const PT Zero(0), One(1);
    const T  COne(1);

    Pair<IndexType> wsQuery = ev_wsq(A);
    IndexType minWork = wsQuery.first;
    IndexType maxWork = wsQuery.second;

//
//  Perform and apply worksize query if requested
//
    if (work.length()!=0 && work.length()<minWork) {
        ASSERT(0);
    } else if (work.length()==0) {
        work.resize(maxWork);
    }

    const IndexType lWork = work.length();
    const IndexType n     = A.dim();

//
//  Quick return if possible
//
    if (n==0) {
        return 0;
    }

    if (n==1) {
        w(1) = real(A(1,1));
        work(1) = 1;
        if (computeV) {
            A(1,1) = COne;
        }
        return 0;
    }
//
//  Get machine constants.
//
    const PT safeMin  = lamch<PT>(SafeMin);
    const PT eps      = lamch<PT>(Precision);
    const PT smallNum = safeMin / eps;
    const PT bigNum   = One / smallNum;
    const PT rMin     = sqrt(smallNum);
    const PT rMax     = sqrt(bigNum);
//
//  Scale matrix to allowable range, if necessary.
//
    const PT ANorm = lan(MaximumNorm, A);
    bool scaleA = false;
    PT sigma;
    if (ANorm>Zero && ANorm<rMin) {
        scaleA = true;
        sigma  = rMin / ANorm;
    } else if (ANorm>rMax) {
        scaleA = true;
        sigma  = rMax / ANorm;
    }
    if (scaleA) {
        lascl(upper ? LASCL::UpperTriangular : LASCL::LowerTriangular,
              IndexType(0), IndexType(0), One, sigma, A);
    }
//
//  Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
//
    auto e_    = rWork(_(1,n-1));
    auto tau_  = work(_(1,n-1));
    auto work_ = work(_(n+1, lWork));
    trd(A, w, e_, tau_, work_);
//
//  For eigenvalues only, call DSTERF.  For eigenvectors, first call
//  ZUNGTR to generate the unitary matrix, then call ZSTEQR.
//
    IndexType info = 0;

    if (!computeV) {
        info = sterf(w, e_);
    } else {
        ungtr(A, tau_, work_);

        auto rWork_ = rWork(_(n+1, 3*n-2));

        STEQR::ComputeZ   jobZ = (computeV)
                               ? STEQR::Orig
                               : STEQR::No;

        steqr(jobZ, w, e_, A.general(), rWork_);
    }
//
//  If matrix was scaled, then rescale eigenvalues appropriately.
//
    if (scaleA) {
        if (info==0) {
            w *= One/sigma;
        } else {
            w(_(1,info-1)) *= One/sigma;
        }
    }
//
//  Set WORK(1) to optimal complex workspace size.
//
    work(1) = maxWork;
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (he)ev [worksize query hermitian variant] ---------------------------------

template <typename MA>
typename RestrictTo<IsComplex<typename MA::ElementType>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq_impl(const HeMatrix<MA>  &A)
{
    using std::max;

    typedef typename HeMatrix<MA>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef typename HeMatrix<MA>::IndexType            IndexType;

//
//  Compute minimal workspace
//
    IndexType  n = A.dim();
    IndexType  minWork;

    if (n==0) {
        minWork = 1;
    } else {
        minWork = max(1,2*n-1);
    }

//
//  Get optimal workspace from external LAPACK
//
    T           DUMMY, WORK;
    PT          RDUMMY, RWORK;
    IndexType   LWORK = -1;
    cxxlapack::heev('N',
                    getF77Char(A.upLo()),
                    A.dim(),
                    &DUMMY,
                    A.leadingDimension(),
                    &RDUMMY,
                    &WORK,
                    LWORK,
                    &RWORK);
    return Pair<IndexType>(minWork,WORK.real());
}

//-- (he)ev [hermitian variant] ------------------------------------------------

template <typename MA, typename VW, typename VWORK, typename VRWORK>
typename HeMatrix<MA>::IndexType
ev_impl(bool                  computeV,
        HeMatrix<MA>          &A,
        DenseVector<VW>       &w,
        DenseVector<VWORK>    &work,
        DenseVector<VRWORK>   &rWork)
{
    using std::max;

    typedef typename HeMatrix<MA>::IndexType  IndexType;

    if (work.length()==0) {
        const auto ws = ev_wsq_impl(A);
        work.resize(ws.second, 1);
    }
    IndexType  info;
    info = cxxlapack::heev(computeV ? 'V' : 'N',
               getF77Char(A.upLo()),
                           A.dim(),
                           A.data(),
                           A.leadingDimension(),
                           w.data(),
                           work.data(),
                           work.length(),
               rWork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (he)ev [complex variant] --------------------------------------------------
template <typename MA, typename VW, typename VWORK, typename VRWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VW>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeV,
   MA       &&A,
   VW       &&w,
   VWORK    &&work,
   VRWORK   &&rWork)
{
    LAPACK_DEBUG_OUT("(he)ev [complex]");

    using std::max;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;

    const IndexType n = A.dim();

#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VW>::Type      VectorW;
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
    typedef typename RemoveRef<VRWORK>::Type  VectorRWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rWork.firstIndex()==1);

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==0 || w.length()==n);

    if (work.length()!=0) {
        ASSERT(work.length()>=max(IndexType(1),2*n-1));
    }
    if (rWork.length()!=0) {
        ASSERT(rWork.length()>=max(IndexType(1),3*n-2));
    }
#   endif

//
//  Resize output arguments if they are empty and needed
//
    if (w.length()==0) {
        w.resize(n, 1);
    }

    if (rWork.length()==0) {
        rWork.resize(max(IndexType(1),3*n-2));
    }
//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView      A_org     = A;
    typename VectorW::NoView      w_org     = w;
    typename VectorWork::NoView   work_org  = work;
    typename VectorRWork::NoView  rWork_org = rWork;
#   endif

//
//  Call implementation
//
    IndexType result = LAPACK_SELECT::ev_impl(computeV, A, w, work, rWork);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView      A_generic     = A;
    typename VectorW::NoView      w_generic     = w;
    typename VectorWork::NoView   work_generic  = work;
    typename VectorRWork::NoView  rWork_generic = rWork;

    A     = A_org;
    w     = w_org;
    work  = work_org;
    rWork = rWork_org;

    IndexType result_ = external::ev_impl(computeV, A, w, work, rWork);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "A_org = " << A_org << std::endl;
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(w_generic, w, " w_generic", "w")) {
        std::cerr << "CXXLAPACK: w_generic = " << w_generic << std::endl;
        std::cerr << "F77LAPACK: w = " << w << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(rWork_generic, rWork, "rWork_generic", "rWork")) {
        std::cerr << "CXXLAPACK: rWork_generic = " << rWork_generic
                  << std::endl;
        std::cerr << "F77LAPACK: rWork = " << rWork << std::endl;
        failed = true;
    }

    if (! isIdentical(result, result_, " result", "result_")) {
        std::cerr << "CXXLAPACK:  result = " << result << std::endl;
        std::cerr << "F77LAPACK: result_ = " << result_ << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "computeV = " << computeV << std::endl;
        ASSERT(0);
    }
#   endif

    return result;
}

//-- (he)ev [worksize query] ---------------------------------------------------

template <typename MA>
typename RestrictTo<IsHeMatrix<MA>::value,
         Pair<typename MA::IndexType> >::Type
ev_wsq(const MA  &A)
{
    LAPACK_DEBUG_OUT("(he)ev_wsq");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
#   endif

//
//  Call implementation
//
    const auto ws = LAPACK_SELECT::ev_wsq_impl(A);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    const auto optWorkSize = external::ev_wsq_impl(A);
    if (! isIdentical(optWorkSize.first, ws.first,
                      "optWorkSize.first", "ws.first"))
    {
        ASSERT(0);
    }
    if (! isIdentical(optWorkSize.second, ws.second,
                      "optWorkSize.second", "ws.second"))
    {
        ASSERT(0);
    }
#   endif

    return ws;
}

//-- (he)ev [real variant with temporary workspace] ----------------------------

template <typename MA, typename VW>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsRealDenseVector<VW>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeV,
   MA       &&A,
   VW       &&w)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rWork;

    return ev(computeV, A, w, work, rWork);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HE_EV_TCC
