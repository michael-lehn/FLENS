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
      SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
     $                  RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HB_EV_H
#define FLENS_LAPACK_HB_EV_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

#ifdef USE_CXXLAPACK

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

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HB_EV_H
