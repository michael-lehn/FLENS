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
      SUBROUTINE DGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM,
     $                   WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK,
     $                   IWORK, LIWORK, BWORK, INFO )
      SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W,
     $                   VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK,
     $                   BWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 */

#ifndef FLENS_LAPACK_IMPL_ESX_H
#define FLENS_LAPACK_IMPL_ESX_H 1

#include <cmath>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)esx ===================================================================
//
//  Real variant
//
template <typename SelectFunction, typename MA, typename IndexType,
          typename VWR, typename VWI, typename MVS,
          typename RCONDE, typename RCONDV,
          typename VWORK, typename VIWORK, typename VBWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsInteger<IndexType>::value
                     && IsRealDenseVector<VWR>::value
                     && IsRealDenseVector<VWI>::value
                     && IsRealGeMatrix<MVS>::value
                     && IsNotComplex<RCONDE>::value
                     && IsNotComplex<RCONDV>::value
                     && IsRealDenseVector<VWORK>::value
                     && IsRealDenseVector<VIWORK>::value
                     && IsRealDenseVector<VBWORK>::value,
             IndexType>::Type
    esx(bool                computeSchurVectors,
        bool                sortEigenvalues,
        SelectFunction      selectFunction,
        SENSE::Sense        sense,
        MA                  &&A,
        IndexType           &sDim,
        VWR                 &&wr,
        VWI                 &&wi,
        MVS                 &&VS,
        RCONDE              &rCondE,
        RCONDV              &rCondV,
        VWORK               &&work,
        VIWORK              &&iWork,
        VBWORK              &&bWork);

#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename SelectFunction, typename MA, typename IndexType,
          typename VW, typename MVS, typename RCONDE, typename RCONDV,
          typename VWORK, typename VRWORK, typename VBWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsInteger<IndexType>::value
                     && IsComplexDenseVector<VW>::value
                     && IsComplexGeMatrix<MVS>::value
                     && IsNotComplex<RCONDE>::value
                     && IsNotComplex<RCONDV>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value
                     && IsRealDenseVector<VBWORK>::value,
             IndexType>::Type
    esx(bool                computeSchurVectors,
        bool                sortEigenvalues,
        SelectFunction      selectFunction,
        SENSE::Sense        sense,
        MA                  &&A,
        IndexType           &sDim,
        VW                  &&w,
        MVS                 &&VS,
        RCONDE              &rCondE,
        RCONDV              &rCondV,
        VWORK               &&work,
        VRWORK              &&rWork,
        VBWORK              &&bWork);

#endif // USE_CXXLAPACK

//== workspace query ===========================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsNotComplex<typename MA::ElementType>::value,
             Quadruple<typename MA::IndexType> >::Type
    esx_wsq(bool                 computeSchurVectors,
            SENSE::Sense         sense,
            const GeMatrix<MA>   &A);

#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplex<typename MA::ElementType>::value,
             Triple<typename MA::IndexType> >::Type
    esx_wsq(bool                 computeSchurVectors,
            SENSE::Sense         sense,
            const GeMatrix<MA>   &A);

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ESX_H
