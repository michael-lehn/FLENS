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
       SUBROUTINE ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_UNGLQ_TCC
#define FLENS_LAPACK_IMPL_UNGLQ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- unglq ---------------------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
unglq_impl(IndexType                k,
           GeMatrix<MA>             &A,
           const DenseVector<VTAU>  &tau,
           DenseVector<VWORK>       &work)
{
    cxxlapack::unglq<IndexType>(A.numRows(),
                                A.numCols(),
                                k,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- unglq ---------------------------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
unglq(IndexType k, MA &&A, const VTAU &tau, VWORK &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   ElementType;
    typedef typename RemoveRef<VWORK>::Type VectorWork;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==IndexType(1));
    ASSERT(A.firstCol()==IndexType(1));
    ASSERT(tau.firstIndex()==IndexType(1));
    ASSERT(tau.length()==k);
    ASSERT((work.length()==0) || (work.length()>=A.numRows()));

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(n>=m);
    ASSERT(m>=k);
    ASSERT(k>=0);
#   endif 

//
//  Call implementation
//
    external::unglq_impl(k, A, tau, work);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNGLQ_TCC
