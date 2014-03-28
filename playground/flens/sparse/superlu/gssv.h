/*
 *   Copyright (c) 2013, Klaus Pototzky
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

#ifndef PLAYGROUND_FLENS_SPARSE_SUPERLU_GSSV_H
#define PLAYGROUND_FLENS_SPARSE_SUPERLU_GSSV_H 1

#ifdef WITH_SUPERLU

namespace flens { namespace superlu {

template <typename T>
    typename RestrictTo<IsSame<float, T>::value, int>::Type
    gssv(superlu_options_t  *options,
         SuperMatrix        *A,
         int                *pc,
         int                *pr,
         SuperMatrix        *L,
         SuperMatrix        *U,
         SuperMatrix        *B,
         SuperLUStat_t      *stat);

template <typename T>
    typename RestrictTo<IsSame<double, T>::value, int>::Type
    gssv(superlu_options_t  *options,
         SuperMatrix        *A,
         int                *pc,
         int                *pr,
         SuperMatrix        *L,
         SuperMatrix        *U,
         SuperMatrix        *B,
         SuperLUStat_t      *stat);

template <typename T>
    typename RestrictTo<IsSame<std::complex<float>, T>::value, int >::Type
    gssv(superlu_options_t  *options,
         SuperMatrix        *A,
         int                *pc,
         int                *pr,
         SuperMatrix        *L,
         SuperMatrix        *U,
         SuperMatrix        *B,
         SuperLUStat_t      *stat);

template <typename T>
    typename RestrictTo<IsSame<std::complex<double>, T>::value, int >::Type
    gssv(superlu_options_t  *options,
         SuperMatrix        *A,
         int                *pc,
         int                *pr,
         SuperMatrix        *L,
         SuperMatrix        *U,
         SuperMatrix        *B,
         SuperLUStat_t      *stat);

} } // namespace superlu, flens

#endif // WITH_SUPERLU

#endif // PLAYGROUND_FLENS_SPARSE_SUPERLU_GSSV_H
