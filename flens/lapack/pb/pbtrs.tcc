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
       SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
       SUBROUTINE ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011
 */

#ifndef FLENS_LAPACK_PB_PBTRS_TCC
#define FLENS_LAPACK_PB_PBTRS_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {


#ifdef USE_CXXLAPACK

namespace external {

//-- pbtrs [real variant] ------------------------------------------------------

template <typename MA, typename MB>
typename SbMatrix<MA>::IndexType
pbtrs_impl(SbMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename SbMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::pbtrs<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                A.numOffDiags(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- pbtrs [complex variant] ---------------------------------------------------

template <typename MA, typename MB>
typename HbMatrix<MA>::IndexType
pbtrs_impl(HbMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename HbMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::pbtrs<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                A.numOffDiags(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                B.data(),
                                                B.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external


//== pbtrs =====================================================================
//
//  Real and complex variant
//
template <typename MA, typename MB>
typename RestrictTo<(IsSbMatrix<MA>::value || IsHbMatrix<MA>::value)
                 && IsGeMatrix<MB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
pbtrs(MA &&A, MB &&B)
{
    LAPACK_DEBUG_OUT("pbtrs [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);

    ASSERT(B.numRows()==A.dim());
#   endif

//
//  Call implementation
//
    const IndexType info = external::pbtrs_impl(A, B);

    return info;
}

//== pbsv variant if rhs is vector =============================================
//
//  Real and complex variant
//
template <typename MA, typename VB>
typename RestrictTo<(IsSbMatrix<MA>::value || IsHbMatrix<MA>::value)
                 && IsDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
pbtrs(MA &&A, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::noViewOrder;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b);

    return pbtrs(A, B);
}

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_PB_PBTRS_TCC
