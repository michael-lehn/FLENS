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
       SUBROUTINE DLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
      $                   WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LAQP2_H
#define FLENS_LAPACK_LA_LAQP2_H 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== laqps =====================================================================
//
//  Real and complex variant
//
template <typename IndexType, typename MA, typename JPIV, typename VTAU,
          typename VN1, typename VN2, typename VWORK>
    typename RestrictTo<((IsRealGeMatrix<MA>::value &&
                          IsRealDenseVector<VTAU>::value &&
                          IsRealDenseVector<VWORK>::value)
                      || (IsComplexGeMatrix<MA>::value &&
                          IsComplexDenseVector<VTAU>::value &&
                          IsComplexDenseVector<VWORK>::value))
                      && IsIntegerDenseVector<JPIV>::value
                      && IsRealDenseVector<VN1>::value
                      && IsRealDenseVector<VN2>::value,
             void>::Type
    laqp2(IndexType  offset,
          MA         &&A,
          JPIV       &&jPiv,
          VTAU       &&tau,
          VN1        &&vn1,
          VN2        &&vn2,
          VWORK      &&work);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQP2_H
