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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_H 1

#include <cxxblas/typedefs.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
    void
    PackMatrixColToColMajor( IndexType k, const T *A, IndexType ldA, T *A_to, IndexType ldA_to);
    
template <typename IndexType, typename T>
    void
    PackMatrixColToColMajor_4( IndexType k, const T *A, IndexType ldA, T *A_to );

template <typename IndexType, typename T>
    void
    PackMatrixColToRowMajor( IndexType k, const T *B, IndexType ldB, T *b_to, IndexType ldB_to );
    
template <typename IndexType, typename T>
    void
    PackMatrixColToRowMajor_4( IndexType k, const T *B, IndexType ldB, T *b_to );
    
template <typename IndexType, typename T>
    void
    PackMatrixColToColMajor_conj( IndexType k, const T *A, IndexType ldA, T *A_to, IndexType ldA_to );
    
template <typename IndexType, typename T>
    void
    PackMatrixColToColMajor_4_conj( IndexType k, const T *A, IndexType ldA, T *A_to );
    
template <typename IndexType, typename T>
    void
    PackMatrixColToRowMajor_conj( IndexType k, const T *B, IndexType ldB, T *b_to, IndexType ldB_to );
    
template <typename IndexType, typename T>
    void
    PackMatrixColToRowMajor_4_conj( IndexType k, const T *B, IndexType ldB, T *b_to );

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_H
