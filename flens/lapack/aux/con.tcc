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

/*
 */

#ifndef FLENS_LAPACK_AUX_LACON_TCC
#define FLENS_LAPACK_AUX_LACON_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
void
con_generic(Norm                norm,
            const GeMatrix<MA>  &A,
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

    const IndexType n = A.numRows();

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
    ElementType normInvA = Zero;
    bool normIn = false;

    IndexType kase, kase1;

    if (norm==OneNorm) {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;

    auto work1 = work(_(1,n));
    auto work2 = work(_(n+1,2*n));
    auto work3 = work(_(2*n+1,3*n));
    auto work4 = work(_(3*n+1,4*n));

    bool computeRCond = true;

    while (true) {
        lacn2(work2, work1, iwork, normInvA, kase, iSave);
        if (kase==0) {
            break;
        }

        ElementType sl, su;
        if (kase==kase1) {
//
//          Multiply by inv(L).
//
            latrs(NoTrans, normIn, A.lowerUnit(), work1, sl, work3);
//
//          Multiply by inv(U).
//
            latrs(NoTrans, normIn, A.upper(), work1, su, work4);
        } else {
//
//          Multiply by inv(U**T).
//
            latrs(Trans, normIn, A.upper(), work1, su, work4);
//
//          Multiply by inv(L**T).
//
            latrs(Trans, normIn, A.lowerUnit(), work1, sl, work3);
        }
//
//      Divide X by 1/(SL*SU) if doing so will not cause overflow.
//
        const ElementType scale = sl*su;
        normIn = true;
        if (scale!=One) {
            const IndexType ix = blas::iamax(work1);
            if (scale<abs(work1(ix))*smallNum || scale==Zero) {
                computeRCond = false;
                break;
            }
            rscl(scale, work1);
        }
    }
//
//  Compute the estimate of the reciprocal condition number.
//
    if (computeRCond) {
        if (normInvA!=Zero) {
            rCond = (One/normInvA) / normA;
        }
    }
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
void
con_native(Norm                norm,
           const GeMatrix<MA>  &A,
           const NORMA         &normA,
           RCOND               &rCond,
           DenseVector<VWORK>  &work,
           DenseVector<VIWORK> &iwork)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const char       NORM  = char(norm);
    const INTEGER    N     = A.numRows();
    const INTEGER    LDA   = A.leadingDimension();
    const T          ANORM = normA;
    INTEGER          INFO;

    if (IsSame<T,double>::value) {
        LAPACK_IMPL(dgecon)(&NORM,
                            &N,
                            A.data(),
                            &LDA,
                            &ANORM,
                            &rCond,
                            work.data(),
                            iwork.data(),
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO>=0);
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
void
con(Norm                norm,
    const GeMatrix<MA>  &A,
    const NORMA         &normA,
    RCOND               &rCond,
    DenseVector<VWORK>  &work,
    DenseVector<VIWORK> &iwork)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(norm==InfinityNorm || norm==OneNorm);
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==4*n);

    ASSERT(iwork.firstIndex()==1);
    ASSERT(iwork.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    RCOND                                rCond_org = rCond;
    typename DenseVector<VWORK>::NoView  work_org  = work;
    typename DenseVector<VIWORK>::NoView iwork_org = iwork;
#   endif

//
//  Call implementation
//
    con_generic(norm, A, normA, rCond, work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    RCOND                                rCond_generic = rCond;
    typename DenseVector<VWORK>::NoView  work_generic  = work;
    typename DenseVector<VIWORK>::NoView iwork_generic = iwork;

    rCond = rCond_org;
    work  = work_org;
    iwork = iwork_org;

    con_native(norm, A, normA, rCond, work, iwork);

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

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename NORMA, typename RCOND,
          typename VWORK, typename VIWORK>
void
con(Norm         norm,
    const MA     &A,
    const NORMA  &normA,
    RCOND        &&rCond,
    VWORK        &&work,
    VIWORK       &&iwork)
{
    CHECKPOINT_ENTER;
    con(norm, A, normA, rCond, work, iwork);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LACON_TCC
