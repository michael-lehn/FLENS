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
       SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
       SUBROUTINE ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_PB_PBTRF_TCC
#define FLENS_LAPACK_PB_PBTRF_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {


#ifdef USE_CXXLAPACK

namespace external {

//-- pbtrf [real variant] ------------------------------------------------------

template <typename MA>
typename SbMatrix<MA>::IndexType
pbtrf_impl(SbMatrix<MA> &A)
{
    typedef typename SbMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::pbtrf<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.numOffDiags(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

//-- pbtrf [complex variant] ---------------------------------------------------

template <typename MA>
typename HbMatrix<MA>::IndexType
pbtrf_impl(HbMatrix<MA> &A)
{
    typedef typename HbMatrix<MA>::IndexType  IndexType;

    IndexType info = cxxlapack::pbtrf<IndexType>(getF77Char(A.upLo()),
                                                 A.dim(),
                                                 A.numOffDiags(),
                                                 A.data(),
                                                 A.leadingDimension());
    ASSERT(info>=0);
    return info;
}

} // namespace external


//== pbtrf =====================================================================
//
//  Real and complex variant
//
template <typename MA>
typename RestrictTo<IsSbMatrix<MA>::value || IsHbMatrix<MA>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
pbtrf(MA &&A)
{
    LAPACK_DEBUG_OUT("pbtrf [complex]");

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
#   endif

//
//  Call implementation
//
    const IndexType info = external::pbtrf_impl(A);

    return info;
}

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_PB_PBTRF_TCC
