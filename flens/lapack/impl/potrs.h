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
       SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
       SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_POTRS_H
#define FLENS_LAPACK_IMPL_POTRS_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== potrs =====================================================================
//
//  Real variant
//
template <typename MA, typename MB>
    typename RestrictTo<IsRealSyMatrix<MA>::value
                     && IsRealGeMatrix<MB>::value,
             void>::Type
    potrs(const MA &A, MB &&B);


#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA, typename MB>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsComplexGeMatrix<MB>::value,
             void>::Type
    potrs(const MA &A, MB &&B);

#endif // USE_CXXLAPACK


//
//  Variant for convenience: Rhs b is vector
//
template <typename MA, typename VB>
    typename RestrictTo<(IsSyMatrix<MA>::value || IsHeMatrix<MA>::value)
                     && IsDenseVector<VB>::value,
             void>::Type
    potrs(const MA &A, VB &&b);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_POTRS_H
