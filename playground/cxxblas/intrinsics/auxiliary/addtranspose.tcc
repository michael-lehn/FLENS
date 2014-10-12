/*
 *   Copyright (c) 2012, Klaus Pototzky
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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_TRANSPOSE_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_TRANSPOSE_TCC 1

#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/classes/classes.h>

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
void
add_transpose_8x8(const T *A, const IndexType ldA , T *B, const IndexType ldB)
{
    for (IndexType i=0; i<8; ++i) {
        for (IndexType j=0; j<8; ++j) {
            B[i*ldB+j] += A[j*ldA+i];
        }
    }
}

#ifdef HAVE_AVX

template <typename IndexType>
void
add_transpose_8x8(const float *A, const IndexType ldA ,
                  float *B, const IndexType ldB)
{
    typedef Intrinsics<float, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    IntrinsicType A_[8], B_[8], tmp_[8];

    for (IndexType k=0; k<8; ++k) {
        A_[k].loadu(A+k*ldA);
    }

    for (IndexType k=0; k<4; ++k) {
        tmp_[2*k  ] = intrinsic_unpacklo_(A_[2*k], A_[2*k+1]);
        tmp_[2*k+1] = intrinsic_unpackhi_(A_[2*k], A_[2*k+1]);
    }

    A_[0] = mm256_shuffle_ps_(tmp_[0].get(),tmp_[2].get(),68);
    A_[1] = mm256_shuffle_ps_(tmp_[0].get(),tmp_[2].get(),238);
    A_[2] = mm256_shuffle_ps_(tmp_[1].get(),tmp_[3].get(),68);
    A_[3] = mm256_shuffle_ps_(tmp_[1].get(),tmp_[3].get(),238);
    A_[4] = mm256_shuffle_ps_(tmp_[4].get(),tmp_[6].get(),68);
    A_[5] = mm256_shuffle_ps_(tmp_[4].get(),tmp_[6].get(),238);
    A_[6] = mm256_shuffle_ps_(tmp_[5].get(),tmp_[7].get(),68);
    A_[7] = mm256_shuffle_ps_(tmp_[5].get(),tmp_[7].get(),238);

    for (IndexType k=0; k<4; ++k) {
        tmp_[k  ] = mm256_permute2f128_ps_(A_[k].get(), A_[k+4].get(), 32);
        tmp_[k+4] = mm256_permute2f128_ps_(A_[k].get(), A_[k+4].get(), 49);
    }


    for (IndexType k=0; k<8; ++k) {
        B_[k].loadu(B+k*ldB);
        B_[k] = intrinsic_add_(B_[k], tmp_[k]);
        B_[k].storeu(B+k*ldB);
    }

}

template <typename IndexType>
void
add_transpose_8x8(const double *A, const IndexType ldA , double *B, const IndexType ldB)
{
    typedef Intrinsics<double, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    IntrinsicType A_[4], B_[4], tmp_[4];

    for (IndexType i=0; i<2; ++i) {
        for (IndexType j=0; j<2; ++j) {

            for (IndexType k=0; k<4; ++k) {
                A_[k].loadu(A+(k+4*i)*ldA+4*j);
            }

            tmp_[0] = mm256_unpacklo_pd_(A_[0].get(), A_[1].get());
            tmp_[1] = mm256_unpackhi_pd_(A_[0].get(), A_[1].get());
            tmp_[2] = mm256_unpacklo_pd_(A_[2].get(), A_[3].get());
            tmp_[3] = mm256_unpackhi_pd_(A_[2].get(), A_[3].get());

            A_[0] = mm256_permute2f128_pd_(tmp_[0].get(), tmp_[2].get(), 32);
            A_[1] = mm256_permute2f128_pd_(tmp_[1].get(), tmp_[3].get(), 32);
            A_[2] = mm256_permute2f128_pd_(tmp_[0].get(), tmp_[2].get(), 49);
            A_[3] = mm256_permute2f128_pd_(tmp_[1].get(), tmp_[3].get(), 49);

            for (IndexType k=0; k<4; ++k) {
                B_[k].loadu(B+(k+4*j)*ldB+4*i);
                B_[k] = intrinsic_add_(B_[k], A_[k]);
                B_[k].storeu(B+(k+4*j)*ldB+4*i);
            }
        }
    }
}

#endif // HAVE_AVX

#endif // USE_INTRINSIC

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_TRANSPOSE_TCC
