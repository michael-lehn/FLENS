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
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
       SUBROUTINE DSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
     $                  INFO )
      SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
     $                  RWORK, INFO )
      SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
      SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,
     $                  INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_EV_H
#define FLENS_LAPACK_IMPL_EV_H 1

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


#ifdef USE_CXXLAPACK

//== (ge)ev ====================================================================
//
//  complex variant
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

//== (he)ev ====================================================================
//
//  complex variant
//
template <typename MA, typename VW, typename VWORK>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsRealDenseVector<VW>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeV,
       MA       &&A,
       VW       &&w,
       VWORK    &&work);

//== (sy)ev ====================================================================
//
//  real variant
//
template <typename MA, typename VW, typename VWORK>
    typename RestrictTo<IsRealSyMatrix<MA>::value
                     && IsRealDenseVector<VW>::value
                     && IsRealDenseVector<VWORK>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeV,
       MA       &&A,
       VW       &&w,
       VWORK    &&work);

//== (hb)ev ====================================================================
//
//  complex variant
//
template <typename MA, typename VW, typename MZ, typename VWORK,
          typename VRWORK>
    typename RestrictTo<IsHbMatrix<MA>::value
                    && IsRealDenseVector<VW>::value
                    && IsComplexGeMatrix<MZ>::value
                    && IsComplexDenseVector<VWORK>::value
                    && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z,
      VWORK    &&work,
      VRWORK   &&rwork);

//== (sb)ev ====================================================================
//
//  Real variant
//
template <typename MA, typename VW, typename MZ, typename VWORK>
    typename RestrictTo<IsRealSbMatrix<MA>::value
                    && IsRealDenseVector<VW>::value
                    && IsRealGeMatrix<MZ>::value
                    && IsRealDenseVector<VWORK>::value,
            typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z,
      VWORK    &&work);

//== (hb)ev ====================================================================
//
//  complex variant
//
template <typename MA, typename VW, typename MZ, typename VWORK,
          typename VRWORK>
    typename RestrictTo<IsHpMatrix<MA>::value
                    && IsRealDenseVector<VW>::value
                    && IsComplexGeMatrix<MZ>::value
                    && IsComplexDenseVector<VWORK>::value
                    && IsRealDenseVector<VRWORK>::value,
            typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z,
      VWORK    &&work,
      VRWORK   &&rwork);

//== (sp)ev ====================================================================
//
//  Real variant
//
template <typename MA, typename VW, typename MZ, typename VWORK>
    typename RestrictTo<IsRealSpMatrix<MA>::value
                    && IsRealDenseVector<VW>::value
                    && IsRealGeMatrix<MZ>::value
                    && IsRealDenseVector<VWORK>::value,
            typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z,
      VWORK    &&work);

#endif // USE_CXXLAPACK


//== (ge)ev workspace query ======================================================
//
//  Real variant
//
template <typename MA>
    typename RestrictTo<IsRealGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, bool computeVR, const MA &A);


#ifdef USE_CXXLAPACK

//== (ge)ev workspace query =======================================================
//
//  complex variant
//
template <typename MA>
    typename RestrictTo<IsComplexGeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, bool computeVR, const MA &A);

//== (he)ev workspace query =======================================================
//
//  complex variant
//
template <typename MA>
    typename RestrictTo<IsHeMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, const MA &A);

//== (sy)ev workspace query =======================================================
//
//  real variant
//
template <typename MA>
    typename RestrictTo<IsRealSyMatrix<MA>::value,
             Pair<typename MA::IndexType> >::Type
    ev_wsq(bool computeVL, const MA &A);

#endif // USE_CXXLAPACK

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


#ifdef USE_CXXLAPACK
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

//== (he)ev ====================================================================
//
//  Complex variant with temporary workspace
//
template <typename MA, typename VW>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsRealDenseVector<VW>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeV,
       MA       &&A,
       VW       &&w);

//== (sy)ev ====================================================================
//
//  real variant with temporary workspace
//
template <typename MA, typename VW>
    typename RestrictTo<IsRealSyMatrix<MA>::value
                     && IsRealDenseVector<VW>::value,
             typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeV,
       MA       &&A,
       VW       &&w);

//== (hb)ev ====================================================================
//
//  complex variant with temporary workspace
//
template <typename MA, typename VW, typename MZ>
    typename RestrictTo<IsHbMatrix<MA>::value
                      && IsRealDenseVector<VW>::value
                      && IsComplexGeMatrix<MZ>::value,
              typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z);


//== (sb)ev ====================================================================
//
//  real variant with temporary workspace
//
template <typename MA, typename VW, typename MZ>
    typename RestrictTo<IsRealSbMatrix<MA>::value
                      && IsRealDenseVector<VW>::value
                      && IsRealGeMatrix<MZ>::value,
              typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z);

//== (hp)ev ====================================================================
//
//  complex variant with temporary workspace
//
template <typename MA, typename VW, typename MZ>
    typename RestrictTo<IsHpMatrix<MA>::value
                      && IsRealDenseVector<VW>::value
                      && IsComplexGeMatrix<MZ>::value,
              typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z);


//== (sp)ev ====================================================================
//
//  real variant with temporary workspace
//
template <typename MA, typename VW, typename MZ>
    typename RestrictTo<IsRealSpMatrix<MA>::value
                      && IsRealDenseVector<VW>::value
                      && IsRealGeMatrix<MZ>::value,
              typename RemoveRef<MA>::Type::IndexType>::Type
    ev(bool     computeZ,
      MA       &&A,
      VW       &&w,
      MZ       &&Z);

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_EV_H
