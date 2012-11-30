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
add_transpose_8x8(const float *A, const IndexType ldA , float *B, const IndexType ldB)
{
    typedef Intrinsics<float, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    IntrinsicType _A[8], _B[8], _tmp[8];

    for (IndexType k=0; k<8; ++k) {
        _A[k].loadu(A+k*ldA);
    }

    for (IndexType k=0; k<4; ++k) {
        _tmp[2*k  ] = _intrinsic_unpacklo(_A[2*k], _A[2*k+1]);
        _tmp[2*k+1] = _intrinsic_unpackhi(_A[2*k], _A[2*k+1]);
    }

    _A[0] = _mm256_shuffle_ps(_tmp[0].get(),_tmp[2].get(),68);
    _A[1] = _mm256_shuffle_ps(_tmp[0].get(),_tmp[2].get(),238);
    _A[2] = _mm256_shuffle_ps(_tmp[1].get(),_tmp[3].get(),68);
    _A[3] = _mm256_shuffle_ps(_tmp[1].get(),_tmp[3].get(),238);
    _A[4] = _mm256_shuffle_ps(_tmp[4].get(),_tmp[6].get(),68);
    _A[5] = _mm256_shuffle_ps(_tmp[4].get(),_tmp[6].get(),238);
    _A[6] = _mm256_shuffle_ps(_tmp[5].get(),_tmp[7].get(),68);
    _A[7] = _mm256_shuffle_ps(_tmp[5].get(),_tmp[7].get(),238);

    for (IndexType k=0; k<4; ++k) {
        _tmp[k  ] = _mm256_permute2f128_ps(_A[k].get(), _A[k+4].get(), 32);
        _tmp[k+4] = _mm256_permute2f128_ps(_A[k].get(), _A[k+4].get(), 49);
    }


    for (IndexType k=0; k<8; ++k) {
        _B[k].loadu(B+k*ldB);
        _B[k] = _intrinsic_add(_B[k], _tmp[k]);
        _B[k].storeu(B+k*ldB);
    }

}

template <typename IndexType>
void
add_transpose_8x8(const double *A, const IndexType ldA , double *B, const IndexType ldB)
{
    typedef Intrinsics<double, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    IntrinsicType _A[4], _B[4], _tmp[4];

    for (IndexType i=0; i<2; ++i) {
        for (IndexType j=0; j<2; ++j) {

            for (IndexType k=0; k<4; ++k) {
                _A[k].loadu(A+(k+4*i)*ldA+4*j);
            }

            _tmp[0] = _mm256_unpacklo_pd(_A[0].get(), _A[1].get());
            _tmp[1] = _mm256_unpackhi_pd(_A[0].get(), _A[1].get());
            _tmp[2] = _mm256_unpacklo_pd(_A[2].get(), _A[3].get());
            _tmp[3] = _mm256_unpackhi_pd(_A[2].get(), _A[3].get());

            _A[0] = _mm256_permute2f128_pd(_tmp[0].get(), _tmp[2].get(), 32);
            _A[1] = _mm256_permute2f128_pd(_tmp[1].get(), _tmp[3].get(), 32);
            _A[2] = _mm256_permute2f128_pd(_tmp[0].get(), _tmp[2].get(), 49);
            _A[3] = _mm256_permute2f128_pd(_tmp[1].get(), _tmp[3].get(), 49);

            for (IndexType k=0; k<4; ++k) {
                _B[k].loadu(B+(k+4*j)*ldB+4*i);
                _B[k] = _intrinsic_add(_B[k], _A[k]);
                _B[k].storeu(B+(k+4*j)*ldB+4*i);
            }
        }
    }
}

#endif // HAVE_AVX

#endif // USE_INTRINSIC

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_AUXILIARY_TRANSPOSE_TCC
