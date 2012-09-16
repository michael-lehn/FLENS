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
      SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_LQ2_TCC
#define FLENS_LAPACK_GE_LQ2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename VTAU, typename VWORK>
void
lq2_impl(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = std::min(m, n);

    const ElementType One(1);

    for (IndexType i=1; i<=k; ++i) {
//
//      Generate elementary reflector H(i) to annihilate A(i,i+1:n)
//
        larfg(n-i+1, A(i,i), A(i,_(std::min(i+1,n),n)), tau(i));

        if (i<m) {
//
//          Apply H(i) to A(i+1:m,i:n) from the right
//
            const ElementType Aii = A(i,i);
            A(i,i) = One;
            auto _work = work(_(1, m-i));
            larf(Right, A(i,_(i,n)), tau(i), A(_(i+1,m), _(i,n)), _work);
            A(i,i) = Aii;
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VTAU, typename VWORK>
void
lq2_impl(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::gelq2<IndexType>(A.numRows(),
                                                 A.numCols(),
                                                 A.data(),
                                                 A.leadingDimension(),
                                                 tau.data(),
                                                 work.data());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename VTAU, typename VWORK>
void
lq2(GeMatrix<MA> &A, DenseVector<VTAU> &tau, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(tau.firstIndex()==1);
    ASSERT(work.firstIndex()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    const IndexType k = std::min(m, n);

    ASSERT(tau.length()>=k);
    ASSERT(work.length()>=m);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename DenseVector<VTAU>::NoView  tau_org    = tau;
    typename DenseVector<VTAU>::NoView  work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lq2_impl(A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Restore output arguments
//
    typename GeMatrix<MA>::NoView       A_generic      = A;
    typename DenseVector<VTAU>::NoView  tau_generic    = tau;
    typename DenseVector<VTAU>::NoView  work_generic   = work;

    A = A_org;
    tau = tau_org;
    work = work_org;
//
//  Compare results
//
    external::lq2_impl(A, tau, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VTAU, typename VWORK>
void
lq2(MA &&A, VTAU &&tau, VWORK &&work)
{
    CHECKPOINT_ENTER;
    lq2(A, tau, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_LQ2_TCC
