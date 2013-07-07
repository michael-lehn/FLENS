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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_H
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/iscomplex.h>
#include <flens/auxiliary/isreal.h>
#include <flens/auxiliary/restrictto.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_n_n(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_c_n(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_t_n(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_ct_n(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);

template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_n_t(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_c_t(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_t_t(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_ct_t(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);
    
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_n_c(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_c_c(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_t_c(IndexType m, IndexType n, IndexType k,
                 const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                 const T &beta,  T *C, IndexType ldC,
                 MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_ct_c(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_n_ct(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_c_ct(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_t_ct(IndexType m, IndexType n, IndexType k,
                  const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                  const T &beta,  T *C, IndexType ldC,
                  MA *A_packed, MB *B_packed, bool first_call);
    
template <typename IndexType, typename T, typename MA, typename MB>
    void
    sub_gemm_ct_ct(IndexType m, IndexType n, IndexType k,
                   const T &alpha, const MA *A, IndexType ldA, const MB *B, IndexType ldB,
                   const T &beta,  T *C, IndexType ldC,
                   MA *A_packed, MB *B_packed, bool first_call);


#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_H
