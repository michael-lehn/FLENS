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
       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
      $                  VS, LDVS, WORK, LWORK, BWORK, INFO )
       SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
      $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_IMPL_ES_H
#define FLENS_LAPACK_IMPL_ES_H 1

#include <cmath>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)es ====================================================================
//
//  Real variant
//
template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename VWORK, typename VBWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsInteger<IndexType>::value
                     && IsRealDenseVector<VWR>::value
                     && IsRealDenseVector<VWI>::value
                     && IsRealGeMatrix<MVS>::value
                     && IsRealDenseVector<VWORK>::value
                     && IsRealDenseVector<VBWORK>::value,
             IndexType>::Type
    es(bool                 computeSchurVectors,
       bool                 sortEigenvalues,
       SelectFunction       select,
       MA                   &&A,
       IndexType            &sDim,
       VWR                  &&wr,
       VWI                  &&wi,
       MVS                  &&VS,
       VWORK                &&work,
       VBWORK               &&bWork);


#ifdef USE_CXXLAPACK

//
//  Complex variant
//
template <typename SelectFunction, typename MA, typename IndexType,
          typename VW, typename MVS, typename VWORK, typename VRWORK,
          typename VBWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsInteger<IndexType>::value
                     && IsComplexDenseVector<VW>::value
                     && IsComplexGeMatrix<MVS>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value
                     && IsRealDenseVector<VBWORK>::value,
             IndexType>::Type
    es(bool                 computeSchurVectors,
       bool                 sortEigenvalues,
       SelectFunction       select,
       MA                   &&A,
       IndexType            &sDim,
       VW                   &&w,
       MVS                  &&VS,
       VWORK                &&work,
       VRWORK               &&rwork,
       VBWORK               &&bWork);

#endif // USE_CXXLAPACK


//== workspace query ===========================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsRealGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    es_wsq(bool     computeSchurVectors,
           const MA &A);


#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    es_wsq(bool     computeSchurVectors,
           const MA &A);

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ES_H
