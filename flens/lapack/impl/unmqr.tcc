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
       SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_UNMQR_TCC
#define FLENS_LAPACK_IMPL_UNMQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- unmqr_wsq [worksize query] ------------------------------------------------

template <typename MA, typename MC>
typename GeMatrix<MC>::IndexType
unmqr_wsq_impl(Side              side,
               Transpose         trans,
               GeMatrix<MA>      &A,
               GeMatrix<MC>      &C)
{
    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    T               WORK, DUMMY;
    const IndexType LWORK   = -1;

    cxxlapack::unmqr<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                &DUMMY,
                                C.data(),
                                C.leadingDimension(),
                                &WORK,
                                LWORK);
    return cxxblas::real(WORK);
}

//-- unmqr ---------------------------------------------------------------------

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
unmqr_impl(Side                       side,
           Transpose                  trans,
           GeMatrix<MA>               &A,
           const DenseVector<VTAU>    &tau,
           GeMatrix<MC>               &C,
           DenseVector<VWORK>         &work)
{
    typedef typename GeMatrix<MC>::IndexType  IndexType;

    if (work.length()==0) {
        IndexType lWork = unmqr_wsq_impl(side, trans, A, C);
        work.resize(lWork);
    }

    cxxlapack::unmqr<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                A.numCols(),
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                C.data(),
                                C.leadingDimension(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- unmqr ---------------------------------------------------------------------

#ifdef USE_CXXLAPACK

template <typename MA, typename VTAU, typename MC, typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexGeMatrix<MC>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
unmqr(Side         side,
      Transpose    trans,
      MA           &&A,
      const VTAU   &tau,
      MC           &&C,
      VWORK        &&work)
{
#   ifndef NDEBUG
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename MatrixC::IndexType     IndexType;

//
//  Test the input parameters
//
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numCols();

    ASSERT(tau.length()==k);

    if (side==Left) {
        ASSERT(A.numRows()==m);
    } else {
        ASSERT(A.numRows()==n);
    }

    if (work.length()>0) {
        if (side==Left) {
            ASSERT(work.length()>=n);
        } else {
            ASSERT(work.length()>=m);
        }
    }
#   endif

//
//  Call implementation
//
    external::unmqr_impl(side, trans, A, tau, C, work);
}

//
//  Variant with temporary workspace
//
template <typename MA, typename VTAU, typename MC>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexGeMatrix<MC>::value,
         void>::Type
unmqr(Side         side,
      Transpose    trans,
      MA           &&A,
      const VTAU   &tau,
      MC           &&C)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    unmqr(side, trans, A, tau, C, work);
}

//
//  Variant for convenience: c is vector
//
template <typename MA, typename VTAU, typename VC, typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VC>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
unmqr(Side         side,
      Transpose    trans,
      MA           &&A,
      const VTAU   &tau,
      VC           &&c,
      VWORK        &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VC>::Type    VectorC;

    typedef typename VectorC::ElementType  ElementType;
    typedef typename VectorC::IndexType    IndexType;

    const IndexType    n     = c.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  C(n, 1, c, n);

    unmqr(side, trans, A, tau, C, work);
}

//
//  Variant for convenience: c is vector and workspace gets created
//                           temporarily.
//
template <typename MA, typename VTAU, typename VC>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexDenseVector<VC>::value,
         void>::Type
unmqr(Side         side,
      Transpose    trans,
      MA           &&A,
      const VTAU   &tau,
      VC           &&c)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VC>::Type    VectorC;

    typedef typename VectorC::ElementType  ElementType;
    typedef typename VectorC::IndexType    IndexType;

    const IndexType    n     = c.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  C(n, 1, c, n);

    unmqr(side, trans, A, tau, C);
}

//-- unmqr_wsq [worksize query] ------------------------------------------------

template <typename MA, typename MC>
typename RestrictTo<IsComplexGeMatrix<MA>::value
                 && IsComplexGeMatrix<MC>::value,
         typename RemoveRef<MC>::Type::IndexType>::Type
unmqr_wsq(Side        side,
          Transpose   trans,
          MA          &&A,
          MC          &&C)
{
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename MatrixC::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numCols();

    if (side==Left) {
        ASSERT(A.numRows()==m);
    } else {
        ASSERT(A.numCols()==n);
    }
#   endif

//
//  Call implementation
//
    const IndexType info = external::unmqr_wsq_impl(side, trans, A, C);

    return info;
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNMQR_TCC
