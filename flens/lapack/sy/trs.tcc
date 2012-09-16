/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_SY_TRS_TCC
#define FLENS_LAPACK_SY_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (sy)trs [real variant] ----------------------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(const SyMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::sytrs<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================


#ifdef USE_CXXLAPACK

//-- (sy)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
        void>::Type
trs(const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(sy)trs [real/complex]");
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    typedef typename RemoveRef<MB>::Type     MatrixB;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

//
//  Call implementation
//
    external::trs_impl(A, piv, B);
}

//-- (sy)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VB>::value,
          void>::Type
trs(const MA &A, const VPIV &piv, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    trs(A, piv, B);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SY_TRS_TCC
