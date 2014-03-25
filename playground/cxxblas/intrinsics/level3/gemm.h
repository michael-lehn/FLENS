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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_H 1

#include <playground/cxxblas/intrinsics/level3/gemm/packmatrix.h>
#include <playground/cxxblas/intrinsics/level3/gemm/subgemm.h>
#include <playground/cxxblas/intrinsics/level3/gemm/kernelgemm.h>

namespace cxxblas {

#ifdef USE_INTRINSIC
    
    
/* Block sizes */
#ifndef BLOCKSIZE_GEMM_M
#   define BLOCKSIZE_GEMM_M 256
#endif
    
#ifndef BLOCKSIZE_GEMM_K
#   define BLOCKSIZE_GEMM_K 128
#endif
    
template <typename IndexType, typename T, typename MA, typename MB>
    typename flens::RestrictTo< flens::IsIntrinsicsCompatible<T>::value &&
                                flens::IsIntrinsicsCompatible<MA>::value &&
                                flens::IsIntrinsicsCompatible<MB>::value &&
                                (flens::IsComplex<T>::value ||
                                    (flens::IsReal<MA>::value && flens::IsReal<MB>::value) ),
                               void>::Type
    gemm(StorageOrder order, Transpose transA, Transpose transB,
         IndexType m, IndexType n, IndexType k,
         const T &alpha,
         const MA *A, IndexType ldA,
         const MB *B, IndexType ldB,
         const T &beta,
         T *C, IndexType ldC);
   
#endif // USE_INTRINSIC
    
}

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_H
