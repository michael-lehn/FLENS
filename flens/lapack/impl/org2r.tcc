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
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_ORG2R_TCC
#define FLENS_LAPACK_IMPL_ORG2R_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
org2r_impl(IndexType                 k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    const T  Zero(0), One(1);
//
//  Quick return if possible
//
    if (n<=0) {
        return;
    }
//
//  Initialise columns k+1:n to columns of the unit matrix
//
    for (IndexType j=k+1; j<=n; ++j) {
        A(_(1,m),j) = Zero;
        A(j,j) = One;
    }

    for (IndexType i=k; i>=1; --i) {
//
//      Apply H(i) to A(i:m,i+1:n) from the left
//
        if (i<n) {
            A(i,i) = One;
            larf(Left, A(_(i,m),i), tau(i), A(_(i,m),_(i+1,n)), work);
        }
        if (i<m) {
            blas::scal(-tau(i), A(_(i+1,m),i));
        }
        A(i,i) = One - tau(i);
//
//      Set A(1:i-1,i) to zero
//
        A(_(1,i-1),i) = Zero;
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
org2r_impl(IndexType                 k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    cxxlapack::org2r<IndexType>(A.numRows(),
                                A.numCols(),
                                k,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VTAU>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
org2r(IndexType         k,
      MA                &&A,
      const VTAU        &tau,
      VWORK             &&work)
{
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename RemoveRef<VWORK>::Type  VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));
    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(work.length()>=A.numCols());

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(n<=m);
    ASSERT(k<=n);
    ASSERT(0<=k);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView     A_    = A;
    typename VectorWork::NoView  work_ = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::org2r_impl(k, A, tau, work);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    external::org2r_impl(k, A_, tau, work_);

    bool failed = false;
    if (! isIdentical(A, A_, " A", "A_")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: A_ = " << A_ << std::endl;
        failed = true;
    }

    if (! isIdentical(work, work_, " work", "work_")) {
        std::cerr << "CXXLAPACK:  work = " << work << std::endl;
        std::cerr << "F77LAPACK: work_ = " << work_ << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: org2r.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: org2r.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ORG2R_TCC
