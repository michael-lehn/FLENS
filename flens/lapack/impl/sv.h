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
       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
      $                  LWORK, INFO )
       SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
      $                  LWORK, INFO )
       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
      $                  LWORK, INFO )

       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

       SUBROUTINE DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
       SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
       SUBROUTINE ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_SV_H
#define FLENS_LAPACK_IMPL_SV_H 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)sv ====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsGeMatrix<MB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);

//== sv variant if rhs is vector ===========================================
//
//  Real and complex
//
template <typename MA, typename VPIV, typename VB>
    typename RestrictTo< (IsGeMatrix<MA>::value ||
                          IsHeMatrix<MA>::value ||
                          IsSyMatrix<MA>::value ||
                          IsGbMatrix<MA>::value ||
                          IsHpMatrix<MA>::value ||
                          IsSpMatrix<MA>::value)
                     && IsIntegerDenseVector<VPIV>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, VB &&b);


#ifdef USE_CXXLAPACK

//== (he)sv ====================================================================
//
//  complex variant
//
template <typename MA, typename VPIV, typename MB, typename VWORK>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexGeMatrix<MB>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B, VWORK && work);

//== (he)sv if rhs is vector ===================================================
//
//  complex variant
//
template <typename MA, typename VPIV, typename VB, typename VWORK>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexDenseVector<VB>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, VB &&b, VWORK && work);

//== (he)sv ====================================================================
//
//  complex variant with temporary workspace
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexGeMatrix<MB>::value ,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);

//== (sy)sv ====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename MB, typename VWORK>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsGeMatrix<MB>::value
                     && IsDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B, VWORK &&work);

//== (sy)sv if rhs is vector ===================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename VB, typename VWORK>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsDenseVector<VB>::value
                     && IsDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, VB &&b, VWORK &&work);

//== (sy)sv ====================================================================
//
//  Real and complex variant with temporary workspace
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsGeMatrix<MB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);


//== (gb)sv ====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsGbMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsGeMatrix<MB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);


//== (hp)sv ====================================================================
//
//  complex variant
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsHpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexGeMatrix<MB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);


//== (sp)sv ====================================================================
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename MB>
    typename RestrictTo<IsSpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsGeMatrix<MB>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    sv(MA &&A, VPIV &&piv, MB &&B);


#endif

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_SV_H
