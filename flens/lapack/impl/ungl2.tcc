/*
 *   Copyright (c) 2011, 2013, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_UNGL2_TCC
#define FLENS_LAPACK_IMPL_UNGL2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
ungl2_impl(IndexType                 k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::ElementType  T;

    const T Zero(0), One(1);

    const Underscore<IndexType> _;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//  Quick return if possible
//
    if (m<=0) {
        return;
    }
//
//  Initialise rows k+1:m to rows of the unit matrix
//
    if (k<m) {
        for (IndexType j=1; j<=n; ++j) {
            A(_(k+1,m),j) = Zero;
            if (j>k && j<=m) {
                A(j,j) = One;
            }
        }
    }

    for (IndexType i=k; i>=1; --i) {
//
//      Apply H(i)**H to A(i:m,i:n) from the left
//
        if (i<n) {
            A(i,_(i+1,n)) = conjugate(A(i,_(i+1,n)));
            if (i<m) {
                A(i,i) = One;
                larf(Right, A(i,_(i,n)), conj(tau(i)), A(_(i+1,m),_(i,n)),
                     work);
            }
            A(i,_(i+1,n)) *= -tau(i);
            A(i,_(i+1,n)) = conjugate(A(i,_(i+1,n)));
        }
        A(i,i) = One - conj(tau(i));
//
//      Set A(i,1:i-1) to zero
//
        A(i,_(1,i-1)) = Zero;
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
ungl2_impl(IndexType                 k,
           GeMatrix<MA>              &A,
           const DenseVector<VTAU>   &tau,
           DenseVector<VWORK>        &work)
{
    IndexType info;
    info = cxxlapack::ungl2<IndexType>(A.numRows(),
                                       A.numCols(),
                                       k,
                                       A.data(),
                                       A.leadingDimension(),
                                       tau.data(),
                                       work.data());
    ASSERT(info==0);
    FAKE_USE_NDEBUG(info);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
ungl2(IndexType                 k,
      GeMatrix<MA>              &A,
      const DenseVector<VTAU>   &tau,
      DenseVector<VWORK>        &work)
{
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(tau.length()==k);
    ASSERT(work.length()>=m);

    ASSERT(n>=m);
    ASSERT(m>=k);
    ASSERT(0<=k);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_      = A;
    typename DenseVector<VWORK>::NoView work_   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::ungl2_impl(k, A, tau, work);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    external::ungl2_impl(k, A_, tau, work_);

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
        std::cerr << "error in: ungl2.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: ungl2.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
ungl2(IndexType k, MA &&A, const VTAU &tau, VWORK &&work)
{
    ungl2(k, A, tau, work);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNGL2_TCC
