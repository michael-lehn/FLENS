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
       SUBROUTINE DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK,
      $                   INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_PO_POCON_TCC
#define FLENS_LAPACK_PO_POCON_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- pocon [real variant] ------------------------------------------------------

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
void
pocon_impl(const SyMatrix<MA>  &A,
           const NORMA         &normA,
           RCOND               &rCond,
           DenseVector<VWORK>  &work,
           DenseVector<VIWORK> &iwork)
{
    using std::abs;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType Zero(0), One(1);

    const Underscore<IndexType>  _;

    const bool       upper = (A.upLo()==Upper);
    const IndexType  n     = A.dim();

//
//  Local Arrays
//
    IndexType iSaveData[3] = {0, 0, 0};
    DenseVectorView<IndexType>
        iSave = typename DenseVectorView<IndexType>::Engine(3, iSaveData, 1);
//
//  Quick return if possible
//
    rCond = Zero;
    if (n==0) {
        rCond = One;
        return;
    } else if (normA==Zero) {
        return;
    }

    const ElementType smallNum = lamch<ElementType>(SafeMin);
//
//  Estimate the norm of inv(A).
//
    ElementType  normInvA = Zero;
    ElementType  scale, scaleL, scaleU;
    IndexType    kase = 0;
    bool         normIn = false;

    auto x     = work(_(1,n));
    auto v     = work(_(n+1,2*n));
    auto _work = work(_(2*n+1,3*n));

    while (true) {
        lacn2(v, x, iwork, normInvA, kase, iSave);

        if (kase==0) {
            break;
        }

        if (upper) {
//
//          Multiply by inv(U**T).
//
            latrs(Trans, normIn, A.triangular(), x, scaleL, _work);
            normIn = true;
//
//          Multiply by inv(U).
//
            latrs(NoTrans, normIn, A.triangular(), x, scaleU, _work);
        } else {
//
//          Multiply by inv(L).
//
            latrs(NoTrans, normIn, A.triangular(), x, scaleL, _work);
            normIn = true;
//
//          Multiply by inv(L**T).
//
            latrs(Trans, normIn, A.triangular(), x, scaleU, _work);
        }
//
//      Multiply by 1/SCALE if doing so will not cause overflow.
//
        scale = scaleL * scaleU;
        if (scale!=One) {
            IndexType ix = blas::iamax(x);
            if (scale<abs(work(ix))*smallNum || scale==Zero) {
                break;
            }
            rscl(scale, x);
        }
    }
//
//  Compute the estimate of the reciprocal condition number.
//
    if (normInvA!=Zero) {
        rCond = (One/normInvA) /normA;
    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- pocon [real variant] ------------------------------------------------------

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VWORK2>
void
pocon_impl(const SyMatrix<MA>  &A,
           const NORMA         &normA,
           RCOND               &rCond,
           DenseVector<VWORK>  &work,
           DenseVector<VWORK2> &work2)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::pocon<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                A.data(),
                                A.leadingDimension(),
                                normA,
                                rCond,
                                work.data(),
                                work2.data());
}

//-- pocon [complex variant] ---------------------------------------------------

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VWORK2>
void
pocon_impl(const HeMatrix<MA>  &A,
           const NORMA         &normA,
           RCOND               &rCond,
           DenseVector<VWORK>  &work,
           DenseVector<VWORK2> &work2)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    cxxlapack::pocon<IndexType>(getF77Char(A.upLo()),
                                A.dim(),
                                A.data(),
                                A.leadingDimension(),
                                normA,
                                rCond,
                                work.data(),
                                work2.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- pocon [real variant] ------------------------------------------------------

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
typename RestrictTo<IsRealSyMatrix<MA>::value
                 && IsNotComplex<NORMA>::value
                 && IsNotComplex<RCOND>::value
                 && IsRealDenseVector<VWORK>::value
                 && IsIntegerDenseVector<VIWORK>::value,
         void>::Type
pocon(const MA      &A,
      const NORMA   &normA,
      RCOND         &rCond,
      VWORK         &&work,
      VIWORK        &&iwork)
{
    typedef typename MA::IndexType            IndexType;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
    typedef typename RemoveRef<VIWORK>::Type  VectorIWork;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==3*n);

    ASSERT(iwork.firstIndex()==1);
    ASSERT(iwork.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    RCOND                           rCond_org = rCond;
    typename VectorWork::NoView     work_org  = work;
    typename VectorIWork::NoView    iwork_org = iwork;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::pocon_impl(A, normA, rCond, work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    RCOND                           rCond_generic = rCond;
    typename VectorWork::NoView     work_generic  = work;
    typename VectorIWork::NoView    iwork_generic = iwork;

    rCond = rCond_org;
    work  = work_org;
    iwork = iwork_org;

    external::pocon_impl(A, normA, rCond, work, iwork);

    bool failed = false;
    if (! isIdentical(rCond_generic, rCond, "rCond_generic", "rCond")) {
        std::cerr << "CXXLAPACK: rCond_generic = "
                  << rCond_generic << std::endl;
        std::cerr << "F77LAPACK: rCond = " << rCond << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = "
                  << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(iwork_generic, iwork, "iwork_generic", "iwork")) {
        std::cerr << "CXXLAPACK: iwork_generic = "
                  << iwork_generic << std::endl;
        std::cerr << "F77LAPACK: iwork = " << iwork << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: con" << std::endl;
    }
#   endif
}

//-- pocon [complex variant] ------------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VRWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsNotComplex<NORMA>::value
                 && IsNotComplex<RCOND>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         void>::Type
pocon(const MA      &A,
      const NORMA   &normA,
      RCOND         &rCond,
      VWORK         &&work,
      VRWORK        &&rwork)
{
    typedef typename MA::IndexType            IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType n = A.dim();

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==2*n);

    ASSERT(rwork.firstIndex()==1);
    ASSERT(rwork.length()==n);
#   endif

//
//  Call implementation
//
    external::pocon_impl(A, normA, rCond, work, rwork);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_PO_POCON_TCC
