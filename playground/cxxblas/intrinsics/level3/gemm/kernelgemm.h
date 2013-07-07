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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/iscomplex.h>
#include <flens/auxiliary/isreal.h>
#include <flens/auxiliary/restrictto.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

// Real kernel
template<typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    kernel_gemm_2numElementsx4(IndexType k,
                               const T &alpha,
                               const T *A, IndexType ldA,
                               const T *B, IndexType ldB,
                               const T &beta,
                               T *C, IndexType ldC );

// Complex kernel
template<typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
    kernel_gemm_2numElementsx4(IndexType k,
                               const T &alpha,
                               const T *A, IndexType ldA,
                               const T *B, IndexType ldB,
                               const T &beta,
                               T *C, IndexType ldC );
    
// Mixed Kernel: A and B real, C complex
template<typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    kernel_gemm_2numElementsx4(IndexType k,
                               const std::complex<T> &alpha,
                               const T *A, IndexType ldA,
                               const T *B, IndexType ldB,
                               const std::complex<T> &beta,
                               std::complex<T> *C, IndexType ldC );
    
// Mixed Kernel: A real, B and C complex
template<typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    kernel_gemm_2numElementsx4(IndexType k,
                               const std::complex<T> &alpha,
                               const T *A, IndexType ldA,
                               const std::complex<T> *B, IndexType ldB,
                               const std::complex<T> &beta,
                               std::complex<T> *C, IndexType ldC );
    
// Mixed Kernel: B real, A and C complex
template<typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    kernel_gemm_2numElementsx4(IndexType k,
                               const std::complex<T> &alpha,
                               const std::complex<T> *A, IndexType ldA,
                               const T *B, IndexType ldB,
                               const std::complex<T> &beta,
                               std::complex<T> *C, IndexType ldC );

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_H
