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
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
      
      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )      
      SUBROUTINE ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      
      SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO )
      SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO )
                  
      SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
      SUBROUTINE ZSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
      SUBROUTINE ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_IMPL_TRI_H
#define FLENS_LAPACK_IMPL_TRI_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)tri ===================================================================
//
//  Real variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);


#ifdef USE_CXXLAPACK

//
//  Complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);

#endif // USE_CXXLAPACK

//
//  Real/complex variant with temporary workspace
//
template <typename MA, typename VPIV>
    typename RestrictTo< (IsGeMatrix<MA>::value ||
                          IsHeMatrix<MA>::value ||
                          IsSyMatrix<MA>::value ||
                          IsHpMatrix<MA>::value ||
                          IsSpMatrix<MA>::value )
                     && IsIntegerDenseVector<VPIV>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv);


//== (tr)tri ===================================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsRealTrMatrix<MA>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA &&A);


#ifdef USE_CXXLAPACK

//== (tr)tri ===================================================================    
//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexTrMatrix<MA>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA &&A);
    
//== (he)tri ===================================================================    
//
//  Complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);
    
//== (sy)tri ===================================================================    
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);
    
//== (tp)tri ===================================================================    
//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexTpMatrix<MA>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA &&A);    
    
//== (hp)tri ===================================================================    
//
//  Complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsHpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);
    
//== (sp)tri ===================================================================    
//
//  Real and complex variant
//
template <typename MA, typename VPIV, typename VWORK>
    typename RestrictTo<IsSpMatrix<MA>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    tri(MA          &&A,
        const VPIV  &piv,
        VWORK       &&work);    
#endif // USE_CXXLAPACK

    

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRI_H
