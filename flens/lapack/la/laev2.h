/*
 *   Copyright (c) 2014, Michael Lehn
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
       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
       SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAEV2_H
#define FLENS_LAPACK_LA_LAEV2_H 1

#include <cxxstd/complex.h>

namespace flens { namespace lapack {

//-- laev2 ---------------------------------------------------------------------
//
//  Real variant
//
template <typename T>
    typename RestrictTo<IsReal<T>::value,
             void>::Type
    laev2(const T &a,
          const T &b,
          const T &c,
          T       &rt1,
          T       &rt2,
          T       &cs1,
          T       &sn1);

//
//  Complex variant
//
template <typename T>
    typename RestrictTo<IsReal<T>::value,
             void>::Type
    laev2(const std::complex<T> &a,
          const std::complex<T> &b,
          const std::complex<T> &c,
          T                     &rt1,
          T                     &rt2,
          T                     &cs1,
          std::complex<T>       &sn1);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAEV2_H
