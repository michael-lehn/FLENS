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
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_GE_EV_H
#define FLENS_LAPACK_GE_EV_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {


//== (ge)ev ====================================================================
//
//  Real variant
//
template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VWR>::value
                     && IsRealDenseVector<VWI>::value
                     && IsRealGeMatrix<MVL>::value
                     && IsRealGeMatrix<MVR>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeVL,
       bool     computeVR,
       MA       &&A,
       VWR      &&wr,
       VWI      &&wi,
       MVL      &&VL,
       MVR      &&VR,
       VWORK    &&work);

//== (ge)ev workspace query ====================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsRealGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, bool computeVR, const MA &A);

//== (ge)ev ====================================================================
//
//  Real variant with temporary workspace
//
template <typename MA, typename VWR, typename VWI, typename MVL, typename MVR,
          typename VWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealDenseVector<VWR>::value
                     && IsRealDenseVector<VWI>::value
                     && IsRealGeMatrix<MVL>::value
                     && IsRealGeMatrix<MVR>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeVL,
       bool     computeVR,
       MA       &&A,
       VWR      &&wr,
       VWI      &&wi,
       MVL      &&VL,
       MVR      &&VR);


//== (ge)ev ====================================================================
//
//  Complex variant
//
template <typename MA, typename VW, typename MVL, typename MVR, typename VWORK,
          typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexDenseVector<VW>::value
                     && IsComplexGeMatrix<MVL>::value
                     && IsComplexGeMatrix<MVR>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeVL,
       bool     computeVR,
       MA       &&A,
       VW       &&w,
       MVL      &&VL,
       MVR      &&VR,
       VWORK    &&work,
       VRWORK   &&rWork);

//== (ge)ev workspace query ====================================================
//
//  Complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, bool computeVR, const MA &A);

//== (ge)ev ====================================================================
//
//  Complex variant with temporary workspace
//
template <typename MA, typename VW, typename MVL, typename MVR>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexDenseVector<VW>::value
                     && IsComplexGeMatrix<MVL>::value
                     && IsComplexGeMatrix<MVR>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeVL,
       bool     computeVR,
       MA       &&A,
       VW       &&w,
       MVL      &&VL,
       MVR      &&VR);

//== (ge)ev ====================================================================
//
//  Real/complex variant that only computes eigenvalues
//
template <typename MA, typename VW>
    typename RestrictTo<(IsComplexGeMatrix<MA>::value
                      && IsComplexDenseVector<VW>::value),
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(MA       &&A,
       VW       &&w);

template <typename MA, typename VWR, typename VWI>
    typename RestrictTo<(IsRealGeMatrix<MA>::value
                      && IsRealDenseVector<VWR>::value
                      && IsRealDenseVector<VWI>::value),
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(MA       &&A,
       VWR      &&wr,
       VWI      &&wi);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_EV_H
