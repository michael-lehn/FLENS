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
      SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
     $                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
     $                   LDU, NV, WV, LDWV, NH, WH, LDWH )
      SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
     $                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
     $                   WV, LDWV, NH, WH, LDWH )
 *
 *  -- LAPACK auxiliary routine (version 3.3.0) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2010
 */

#ifndef FLENS_LAPACK_LA_LAQR5_H
#define FLENS_LAPACK_LA_LAQR5_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== laqr5 =====================================================================
//
//  Real variant
//
template <typename IndexType, typename VSR, typename VSI, typename MH,
          typename MZ, typename MV, typename MU, typename MWV, typename MWH>
    typename RestrictTo<IsRealDenseVector<VSR>::value
                     && IsRealDenseVector<VSI>::value
                     && IsRealGeMatrix<MH>::value
                     && IsRealGeMatrix<MZ>::value
                     && IsRealGeMatrix<MV>::value
                     && IsRealGeMatrix<MU>::value
                     && IsRealGeMatrix<MWV>::value
                     && IsRealGeMatrix<MWH>::value,
             void>::Type
    laqr5(bool                      wantT,
          bool                      wantZ,
          IndexType                 kacc22,
          IndexType                 kTop,
          IndexType                 kBot,
          IndexType                 nShifts,
          VSR                       &&sr,
          VSI                       &&si,
          MH                        &&H,
          IndexType                 iLoZ,
          IndexType                 iHiZ,
          MZ                        &&Z,
          MV                        &&V,
          MU                        &&U,
          MWV                       &&WV,
          MWH                       &&WH);

//
//  Complex variant
//
template <typename IndexType, typename VS, typename MH, typename MZ,
          typename MV, typename MU, typename MWV, typename MWH>
    typename RestrictTo<IsComplexDenseVector<VS>::value
                     && IsComplexGeMatrix<MH>::value
                     && IsComplexGeMatrix<MZ>::value
                     && IsComplexGeMatrix<MV>::value
                     && IsComplexGeMatrix<MU>::value
                     && IsComplexGeMatrix<MWV>::value
                     && IsComplexGeMatrix<MWH>::value,
             void>::Type
    laqr5(bool                      wantT,
          bool                      wantZ,
          IndexType                 kacc22,
          IndexType                 kTop,
          IndexType                 kBot,
          IndexType                 nShifts,
          VS                        &&s,
          MH                        &&H,
          IndexType                 iLoZ,
          IndexType                 iHiZ,
          MZ                        &&Z,
          MV                        &&V,
          MU                        &&U,
          MWV                       &&WV,
          MWH                       &&WH);


} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQR5_H
