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
       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )

       SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
       SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
       SUBROUTINE ZHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
       SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )

       SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
       SUBROUTINE ZSPTRF( UPLO, N, AP, IPIV, INFO )
       SUBROUTINE ZHPTRF( UPLO, N, AP, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_TRF_H
#define FLENS_LAPACK_IMPL_TRF_H 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)trf =====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv);


#ifdef USE_CXXLAPACK

//== (he)trf =====================================================================
//
//  complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv, VWORK &&work);

//== (sy)trf =====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv, VWORK &&work);

//== (he/sy)trf with temporary workspace =========================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV>
    typename RestrictTo<(IsHeMatrix<MA>::value ||
                         IsSyMatrix<MA>::value )
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv);

//== (gb)trf =====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV>
    typename RestrictTo<IsGbMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv);


//== (hp)trf =====================================================================
//
// complex variant
//
template <typename MA, typename VPIV>
    typename RestrictTo<IsHpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv);

//== (sp)trf =====================================================================
//
// real and complex variant
//
template <typename MA, typename VPIV>
    typename RestrictTo<IsSpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    trf(MA &&A, VPIV &&piv);

#endif

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRF_H
