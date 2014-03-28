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

#ifndef PLAYGROUND_FLENS_SPARSE_MKL_SV_H
#define PLAYGROUND_FLENS_SPARSE_MKL_SV_H 1

#ifdef WITH_MKLDSS

#include<flens/auxiliary/auxiliary.h>
#include<flens/matrixtypes/matrixtypes.h>
#include<flens/vectortypes/vectortypes.h>

namespace flens { namespace mkldss {

// Interface to MKL DSS, solves AX = B

// General Sparse Matrix
template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsRealGeCRSMatrix<MA>::value &&
               IsRealGeMatrix<MX>::value &&
               IsRealGeMatrix<MB>::value,
               void>::Type
    sv(Transpose transA,
       const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsRealGeCCSMatrix<MA>::value &&
               IsRealGeMatrix<MX>::value &&
               IsRealGeMatrix<MB>::value,
               void>::Type
    sv(Transpose transA,
       const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsComplexGeCRSMatrix<MA>::value &&
              IsComplexGeMatrix<MX>::value &&
              IsComplexGeMatrix<MB>::value,
              void>::Type
    sv(Transpose transA,
       const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsComplexGeCCSMatrix<MA>::value &&
               IsComplexGeMatrix<MX>::value &&
               IsComplexGeMatrix<MB>::value,
               void>::Type
    sv(Transpose transA,
       const MA  &A,
       MX        &&X,
       const MB  &B);

// Symmetric Sparse Matrix
template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsRealSyCCSMatrix<MA>::value &&
               IsRealGeMatrix<MX>::value &&
               IsRealGeMatrix<MB>::value,
               void>::Type
    sv(const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsRealSyCRSMatrix<MA>::value &&
               IsRealGeMatrix<MX>::value &&
               IsRealGeMatrix<MB>::value,
               void>::Type
    sv(const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsComplexSyCCSMatrix<MA>::value &&
               IsComplexGeMatrix<MX>::value &&
               IsComplexGeMatrix<MB>::value,
               void>::Type
    sv(const MA  &A,
       MX        &&X,
       const MB  &B);

template <typename MA, typename MX, typename MB>
    typename
    RestrictTo<IsComplexSyCRSMatrix<MA>::value &&
               IsComplexGeMatrix<MX>::value &&
               IsComplexGeMatrix<MB>::value,
               void>::Type
    sv(const MA  &A,
       MX        &&X,
       const MB  &B);


// Interface for vectors
template <typename MA, typename VX, typename VB>
    typename
    RestrictTo<(IsGeCCSMatrix<MA>::value || IsGeCRSMatrix<MA>::value) &&
               IsDenseVector<VX>::value &&
               IsDenseVector<VB>::value,
               void>::Type
    sv(Transpose transA,
       const MA  &A,
       VX        &&x,
       const VB  &b);

template <typename MA, typename VX, typename VB>
    typename
    RestrictTo<(IsSyCCSMatrix<MA>::value || IsSyCRSMatrix<MA>::value) &&
               IsDenseVector<VX>::value &&
               IsDenseVector<VB>::value,
               void>::Type
    sv(const MA  &A,
       VX        &&x,
       const VB  &b);

} } // namespace mkldss, flens

#endif // WITH_MKLDSS

#endif // PLAYGROUND_FLENS_SPARSE_MKL_SV_H
