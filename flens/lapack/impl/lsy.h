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
       SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_LSY_H
#define FLENS_LAPACK_IMPL_LSY_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)lsy ===================================================================
//
//  Real variant
//
template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealGeMatrix<MB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    lsy(MA           &&A,
        MB           &&B,
        VJPIV        &&jPiv,
        RCOND        rCond,
        VWORK        &&work);



#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK, typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexGeMatrix<MB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             void>::Type
    lsy(MA           &&A,
        MB           &&B,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank,
        VWORK        &&work,
        VRWORK       &&rwork);

#endif // USE_CXXLAPACK


//
//  Real variant with temporary workspace
//
template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealGeMatrix<MB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value,
             void>::Type
    lsy(MA           &&A,
        MB           &&B,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank);


#ifdef USE_CXXLAPACK
//
//  Complex variant with temporary workspace
//
template <typename MA, typename MB, typename VJPIV, typename RCOND,
          typename RANK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexGeMatrix<MB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value,
             void>::Type
    lsy(MA           &&A,
        MB           &&B,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank);

#endif // USE_CXXLAPACK


//== (ge)lsy variant if rhs is vector ==========================================
//
//  Real variant
//
template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value
                     && IsRealDenseVector<VWORK>::value,
             void>::Type
    lsy(MA           &&A,
        VB           &&b,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank,
        VWORK        &&work);


#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename RANK, typename VWORK, typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexDenseVector<VB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             void>::Type
    lsy(MA           &&A,
        VB           &&b,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank,
        VWORK        &&work,
        VRWORK       &&rwork);

#endif // USE_CXXLAPACK

//
//  Real variant with temporary workspace
//
template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename RANK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value,
             void>::Type
    lsy(MA           &&A,
        VB           &&b,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank);


#ifdef USE_CXXLAPACK
//
//  Complex variant with temporary workspace
//
template <typename MA, typename VB, typename VJPIV, typename RCOND,
          typename RANK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexDenseVector<VB>::value
                     && IsIntegerDenseVector<VJPIV>::value
                     && IsReal<RCOND>::value
                     && IsInteger<RANK>::value,
             void>::Type
    lsy(MA           &&A,
        VB           &&b,
        VJPIV        &&jPiv,
        RCOND        rCond,
        RANK         &rank);

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_LSY_H
