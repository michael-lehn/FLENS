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
       SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
      $                  INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_LS_TCC
#define FLENS_LAPACK_IMPL_LS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename MB, typename VWORK>
typename GeMatrix<MA>::IndexType
ls_impl(Transpose                 trans,
        GeMatrix<MA>              &A,
        GeMatrix<MB>              &B,
        DenseVector<VWORK>        &work)
{
    using std::max;
    using flens::min;

    typedef typename GeMatrix<MB>::ElementType  ElementType;
    typedef typename GeMatrix<MB>::IndexType    IndexType;

    const ElementType  Zero(0), One(1);

    const Underscore<IndexType>  _;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  nRhs = B.numCols();
    const IndexType  mn = min(m,n);

    IndexType info = 0;
//
//  Figure out optimal block size
//
    IndexType  nb;
    bool       tpsd = (trans==NoTrans) ? false : true;

    if (m>=n) {
        nb = ilaenv<ElementType>(1, "GEQRF", "", m, n);
        if (tpsd) {
            nb = max(nb, ilaenv<ElementType>(1, "ORMQR", "LN", m, nRhs, n));
        } else {
            nb = max(nb, ilaenv<ElementType>(1, "ORMQR", "LT", m, nRhs, n));
        }
    } else {
        nb = ilaenv<ElementType>(1, "GELQF", "", m, n);
        if (tpsd) {
            nb = max(nb, ilaenv<ElementType>(1, "ORMLQ", "LT", n, nRhs, m));
        } else {
            nb = max(nb, ilaenv<ElementType>(1, "ORMLQ", "LN", n, nRhs, m));
        }
    }

    IndexType  wSize = max(IndexType(1), mn+max(mn,nRhs)*nb);
    if (work.length()==0) {
        work.resize(wSize);
    }
    work(1) = ElementType(wSize);
    const IndexType  lWork = work.length();
//
//  Quick return if possible
//
    if (min(m, n, nRhs)==IndexType(0)) {
        B = Zero;
        return info;
    }
//
//  Get machine parameters
//
    ElementType  smallNum = lamch<ElementType>(SafeMin)
                          / lamch<ElementType>(Precision);
    ElementType  bigNum = One / smallNum;
    labad(smallNum, bigNum);
//
//  Scale A, B if max element outside range [SMLNUM,BIGNUM]
//
    const ElementType  normA = lan(MaximumNorm, A);
    IndexType iScaleA = 0;

    if (normA>Zero && normA<smallNum) {
//
//      Scale matrix norm up to SMLNUM
//
        lascl(LASCL::FullMatrix, 0, 0, normA, smallNum, A);
        iScaleA = 1;
    } else if (normA>bigNum) {
//
//      Scale matrix norm down to BIGNUM
//
        lascl(LASCL::FullMatrix, 0, 0, normA, bigNum, A);
        iScaleA = 2;
    } else if (normA==Zero) {
//
//      Matrix all zero. Return zero solution.
//
        lascl(LASCL::FullMatrix, max(m,n), nRhs, Zero, Zero, A);
        work(1) = wSize;
    }

    auto _B = (tpsd) ? B(_(1,n),_) : B(_(1,m),_);
    const ElementType  normB = lan(MaximumNorm, _B);
    IndexType iScaleB = 0;
    IndexType scaleLen = 0;

    if (normB>Zero && normB<smallNum) {
//
//      Scale matrix norm up to SMLNUM
//
        lascl(LASCL::FullMatrix, 0, 0, normB, smallNum, B);
        iScaleB = 1;
    } else if (normB>bigNum) {
//
//      Scale matrix norm down to BIGNUM
//
        lascl(LASCL::FullMatrix, 0, 0, normB, bigNum, B);
        iScaleB = 2;
    }

    auto tau   = work(_(1,mn));
    auto _work = work(_(mn+1,lWork));
    if (m>=n) {
//
//      compute QR factorization of A
//
        qrf(A, tau, _work);
        const auto R = A(_(1,n),_(1,n)).upper();
//
//      workspace at least N, optimally N*NB
//
        if (!tpsd) {
//
//          Least-Squares Problem min || A * X - B ||
//
//          B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
//
            ormqr(Left, Trans, A, tau, B, _work);
//
//          workspace at least NRHS, optimally NRHS*NB
//
//          B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
//
            info = trs(NoTrans, R, B(_(1,n),_));

            if (info>0) {
                return info;
            }
            scaleLen = n;

        } else {
//
//          Overdetermined system of equations A**T * X = B
//
//          B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
//
            info = trs(Trans, R, B(_(1,n),_));

            if (info>0) {
                return info;
            }
//
//          B(N+1:M,1:NRHS) = ZERO
//
            B(_(n+1,m),_) = Zero;
//
//          B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
//
            ormqr(Left, NoTrans, A, tau, B, _work);
//
//          workspace at least NRHS, optimally NRHS*NB
//
            scaleLen = m;

        }

    } else {
//
//      Compute LQ factorization of A
//
        lqf(A, tau, _work);
        const auto L = A(_(1,m),_(1,m)).lower();
//
//       workspace at least M, optimally M*NB.
//
        if (!tpsd) {
//
//          underdetermined system of equations A * X = B
//
//          B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
//
            info = trs(NoTrans, L, B(_(1,m),_));

            if (info>0) {
                return info;
            }
//
//          B(M+1:N,1:NRHS) = 0
//
            B(_(m+1,n),_) = Zero;
//
//          B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
//
            ormlq(Left, Trans, A, tau, B, _work);
//
//          workspace at least NRHS, optimally NRHS*NB
//
            scaleLen = n;
//
        } else {
//
//          overdetermined system min || A**T * X - B ||
//
//          B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
//
            ormlq(Left, NoTrans, A, tau, B, _work);
//
//          workspace at least NRHS, optimally NRHS*NB
//
//          B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
//
            info = trs(Trans, L, B(_(1,m),_));

            if (info>0) {
                return info;
            }

            scaleLen = m;

        }

    }
//
//  Undo scaling
//
    auto __B = B(_(1,scaleLen),_);
    if (iScaleA==1) {
        lascl(LASCL::FullMatrix, 0, 0, normA, smallNum, __B);
    } else if (iScaleA==2) {
        lascl(LASCL::FullMatrix, 0, 0, normA, bigNum, __B);
    }
    if (iScaleB==1) {
        lascl(LASCL::FullMatrix, 0, 0, smallNum, normB, __B);
    } else if (iScaleB==2) {
        lascl(LASCL::FullMatrix, 0, 0, bigNum, normB, __B);
    }

    work(1) = ElementType(wSize);
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename MB, typename VWORK>
typename GeMatrix<MA>::IndexType
ls_impl(Transpose                 trans,
        GeMatrix<MA>              &A,
        GeMatrix<MB>              &B,
        DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gels<IndexType>(getF77Char(trans),
                                                A.numRows(),
                                                A.numCols(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                B.data(),
                                                B.leadingDimension(),
                                                work.data(),
                                                work.length());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename MB, typename VWORK>
typename GeMatrix<MA>::IndexType
ls(Transpose                 trans,
   GeMatrix<MA>              &A,
   GeMatrix<MB>              &B,
   DenseVector<VWORK>        &work)
{
    using std::max;
    using std::min;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType nRhs = B.numCols();

    ASSERT(B.numRows()==max(m,n));

    if (work.length()>0) {
        const IndexType mn = min(m, n);
        ASSERT(work.length()>=max(IndexType(1),mn+max(mn,nRhs)));
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename GeMatrix<MB>::NoView       B_org      = B;
    typename DenseVector<VWORK>::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::ls_impl(trans, A, B, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic       = A;
    typename GeMatrix<MB>::NoView       B_generic       = B;
    typename DenseVector<VWORK>::NoView work_generic    = work;

//
//  restore output arguments
//
    A = A_org;
    B = B_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::ls_impl(trans, A, B, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
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
        std::cerr << "error in: ls.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: ls.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename MB, typename VWORK>
typename MA::IndexType
ls(Transpose               trans,
   MA                      &&A,
   MB                      &&B,
   VWORK                   &&work)
{
    typename MA::IndexType info;

    CHECKPOINT_ENTER;
    info = ls(trans, A, B, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LS_TCC
