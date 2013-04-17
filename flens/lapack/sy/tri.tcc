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

/* Baesed on
 *
      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_SY_TRI_TCC
#define FLENS_LAPACK_SY_TRI_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (sy)tri [real variant] ----------------------------------------

template <typename MA, typename VP, typename VWORK>
typename SyMatrix<MA>::IndexType
tri_impl(SyMatrix<MA>             &A,
         const DenseVector<VP>    &piv,
         DenseVector<VWORK>       &work)
{
    typedef typename SyMatrix<MA>::IndexType    IndexType;

    IndexType info;

    if (work.length()==0) {
        work.resize(A.dim());
    }

    info = cxxlapack::sytri<IndexType>(getF77Char(A.upLo()),
                                       A.dim(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       work.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

#ifdef USE_CXXLAPACK

//-- (sy)tri -------------------------------------------------------------------

template <typename MA, typename VPIV, typename VWORK>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
tri(MA &&A, const VPIV &piv, VWORK &&work)
{
    using std::max;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

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

    const bool lQuery = (work.length()==0);
    ASSERT(lQuery || work.length()>=max(IndexType(1),n));
#   endif

//
//  Call implementation
//
    const IndexType info = external::tri_impl(A, piv, work);

    return info;
}

//-- (sy)tri [real/complex variant with temporary workspace] -------------------
template <typename MA, typename VPIV>
typename RestrictTo<IsSyMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
tri(MA          &&A,
    const VPIV  &piv)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;
    return tri(A, piv, work);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_SY_TRI_TCC
