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
      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
     $                   VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
      SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     $                   LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE,
     $                   RCONDV, WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_EVX_H
#define FLENS_LAPACK_IMPL_EVX_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)evx ===================================================================
//
//  Real variant
//
template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename VRCONDE, typename VRCONDV, typename VWORK, typename VIWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VWR>::value
                     && IsRealDenseVector<VWI>::value
                     && IsRealGeMatrix<MVL>::value
                     && IsRealGeMatrix<MVR>::value
                     && IsInteger<IndexType>::value
                     && IsRealDenseVector<VSCALE>::value
                     && IsNotComplex<ABNORM>::value
                     && IsRealDenseVector<VRCONDE>::value
                     && IsRealDenseVector<VRCONDV>::value
                     && IsRealDenseVector<VWORK>::value
                     && IsIntegerDenseVector<VIWORK>::value,
             IndexType>::Type
    evx(BALANCE::Balance     balance,
        bool                 computeVL,
        bool                 computeVR,
        SENSE::Sense         sense,
        MA                   &&A,
        VWR                  &&wr,
        VWI                  &&wi,
        MVL                  &&VL,
        MVR                  &&VR,
        IndexType            &iLo,
        IndexType            &iHi,
        VSCALE               &&scale,
        ABNORM               &abNorm,
        VRCONDE              &&rCondE,
        VRCONDV              &&rCondV,
        VWORK                &&work,
        VIWORK               &&iWork);


#ifdef USE_CXXLAPACK

//
//  Complex variant
//
template <typename MA, typename VW, typename MVL, typename MVR,
          typename IndexType, typename VSCALE, typename ABNORM,
          typename VRCONDE, typename VRCONDV, typename VWORK, typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexDenseVector<VW>::value
                     && IsComplexGeMatrix<MVL>::value
                     && IsComplexGeMatrix<MVR>::value
                     && IsInteger<IndexType>::value
                     && IsRealDenseVector<VSCALE>::value
                     && IsNotComplex<ABNORM>::value
                     && IsRealDenseVector<VRCONDE>::value
                     && IsRealDenseVector<VRCONDV>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             IndexType>::Type
    evx(BALANCE::Balance     balance,
        bool                 computeVL,
        bool                 computeVR,
        SENSE::Sense         sense,
        MA                   &&A,
        VW                   &&w,
        MVL                  &&VL,
        MVR                  &&VR,
        IndexType            &iLo,
        IndexType            &iHi,
        VSCALE               &&scale,
        ABNORM               &abNorm,
        VRCONDE              &&rCondE,
        VRCONDV              &&rCondV,
        VWORK                &&work,
        VRWORK               &&rWork);

#endif // USE_CXXLAPACK


//== workspace query ===========================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsRealGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    evx_wsq(bool computeVL, bool computeVR, SENSE::Sense sense, const MA &A);


#ifdef USE_CXXLAPACK

//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    evx_wsq(bool computeVL, bool computeVR, SENSE::Sense sense, const MA &A);

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_EVX_H
