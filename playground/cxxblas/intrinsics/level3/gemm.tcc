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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType>
void
gemm_real_n_n_block_8x8(const float &alpha,
                         const float *A, IndexType ldA,
                         const float *B, IndexType ldB,
                         const float &beta,
                         float *C, IndexType ldC)
{

    CXXBLAS_DEBUG_OUT("gemm_real_n_n_block_8x8 [float]");

    typedef Intrinsics<float, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    ASSERT( beta==float(1) );

    if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::AVX ) {

        IntrinsicType _A;
        IntrinsicType _B[8];
        IntrinsicType _C[8];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<8; ++i) {
            _B[i].loadu(B+i*ldB);
            _B[i] = _intrinsic_mul(_alpha, _B[i]);
        }

        for (IndexType i=0; i<8; ++i) {

            _C[i].loadu(C+i*ldC);

            for (IndexType j=0; j<8; ++j) {
                _A.fill(A[i*ldA+j]);
                _C[i] = _intrinsic_add(_C[i], _intrinsic_mul(_A, _B[j]));
            }

            _C[i].storeu(C+i*ldC);
        }

    } else if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::SSE ) {

        IntrinsicType _A, _C0, _C1;
        IntrinsicType _B[16];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<8; ++i) {
            _B[2*i  ].loadu(B+i*ldB);
            _B[2*i+1].loadu(B+4+i*ldB);

            _B[2*i  ] = _intrinsic_mul(_alpha, _B[2*i  ]);
            _B[2*i+1] = _intrinsic_mul(_alpha, _B[2*i+1]);
        }

        for (IndexType i=0; i<8; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+4);

            for (IndexType j=0; j<8; ++j) {

                _A.fill(A[i*ldA+j]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[2*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[2*j+1]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+4);
        }
    }
}

template <typename IndexType>
void
gemm_real_n_n_block_8x8(const double &alpha,
                         const double *A, IndexType ldA,
                         const double *B, IndexType ldB,
                         const double &beta,
                         double *C, IndexType ldC)
{

    CXXBLAS_DEBUG_OUT("gemm_real_n_n_block_8x8 [double]");

    typedef Intrinsics<double, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;

    ASSERT( beta==double(1) );

    if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::AVX ) {

        IntrinsicType _A, _C0, _C1;
        IntrinsicType _B[16];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<8; ++i) {
            _B[2*i  ].loadu(B+i*ldB);
            _B[2*i+1].loadu(B+4+i*ldB);

            _B[2*i  ] = _intrinsic_mul(_alpha, _B[2*i  ]);
            _B[2*i+1] = _intrinsic_mul(_alpha, _B[2*i+1]);
        }

        for (IndexType i=0; i<8; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+4);

            for (IndexType j=0; j<8; ++j) {

                _A.fill(A[i*ldA+j]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[2*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[2*j+1]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+4);
        }

    } else if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::SSE ) {

        IntrinsicType _A, _C0, _C1, _C2, _C3;
        IntrinsicType _B[32];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<8; ++i) {
            _B[4*i  ].loadu(B+i*ldB);
            _B[4*i+1].loadu(B+i*ldB+2);
            _B[4*i+2].loadu(B+i*ldB+4);
            _B[4*i+3].loadu(B+i*ldB+6);

            _B[4*i  ] = _intrinsic_mul(_alpha, _B[4*i  ]);
            _B[4*i+1] = _intrinsic_mul(_alpha, _B[4*i+1]);
            _B[4*i+2] = _intrinsic_mul(_alpha, _B[4*i+2]);
            _B[4*i+3] = _intrinsic_mul(_alpha, _B[4*i+3]);
        }

        for (IndexType i=0; i<8; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+2);
            _C2.loadu(C+i*ldC+4);
            _C3.loadu(C+i*ldC+6);

            for (IndexType j=0; j<8; ++j) {

                _A.fill(A[i*ldA+j]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[4*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[4*j+1]));
                _C2 = _intrinsic_add(_C2, _intrinsic_mul(_A, _B[4*j+2]));
                _C3 = _intrinsic_add(_C3, _intrinsic_mul(_A, _B[4*j+3]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+2);
            _C2.storeu(C+i*ldC+4);
            _C3.storeu(C+i*ldC+6);
        }
    }
}


template <typename IndexType>
void
gemm_real_t_n_block_8x8(const float &alpha,
                         const float *A, IndexType ldA,
                         const float *B, IndexType ldB,
                         const float &beta,
                         float *C, IndexType ldC)
{

    CXXBLAS_DEBUG_OUT("gemm_real_t_n_block_8x8 [float]");

    typedef Intrinsics<float, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    ASSERT( beta==float(1) );

    if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::AVX ) {

        IntrinsicType _A;
        IntrinsicType _B[numElements];
        IntrinsicType _C;

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<numElements; ++i) {
            _B[i].loadu(B+i*ldB);
            _B[i] = _intrinsic_mul(_alpha, _B[i]);
        }

        for (IndexType i=0; i<numElements; ++i) {

            _C.loadu(C+i*ldC);

            for (IndexType j=0; j<numElements; ++j) {
                _A.fill(A[j*ldA+i]);
                _C = _intrinsic_add(_C, _intrinsic_mul(_A, _B[j]));
            }

            _C.storeu(C+i*ldC);
        }

    } else if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::SSE ) {

        IntrinsicType _A, _C0, _C1;
        IntrinsicType _B[4*numElements];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<2*numElements; ++i) {
            _B[2*i  ].loadu(B+i*ldB);
            _B[2*i+1].loadu(B+numElements+i*ldB);

            _B[2*i  ] = _intrinsic_mul(_alpha, _B[2*i  ]);
            _B[2*i+1] = _intrinsic_mul(_alpha, _B[2*i+1]);
        }

        for (IndexType i=0; i<2*numElements; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+numElements);

            for (IndexType j=0; j<2*numElements; ++j) {

                _A.fill(A[j*ldA+i]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[2*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[2*j+1]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+numElements);
        }
    }
}

template <typename IndexType>
void
gemm_real_t_n_block_8x8(const double &alpha,
                         const double *A, IndexType ldA,
                         const double *B, IndexType ldB,
                         const double &beta,
                         double *C, IndexType ldC)
{

    CXXBLAS_DEBUG_OUT("gemm_real_t_n_block_8x8 [double]");

    typedef Intrinsics<double, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    ASSERT( beta==double(1) );

    if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::AVX ) {

        IntrinsicType _A, _C0, _C1;
        IntrinsicType _B[4*numElements];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<2*numElements; ++i) {
            _B[2*i  ].loadu(B+i*ldB);
            _B[2*i+1].loadu(B+numElements+i*ldB);

            _B[2*i  ] = _intrinsic_mul(_alpha, _B[2*i  ]);
            _B[2*i+1] = _intrinsic_mul(_alpha, _B[2*i+1]);
        }

        for (IndexType i=0; i<2*numElements; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+numElements);

            for (IndexType j=0; j<2*numElements; ++j) {

                _A.fill(A[j*ldA+i]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[2*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[2*j+1]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+numElements);
        }

    } else if ( DEFAULT_INTRINSIC_LEVEL==IntrinsicsLevel::SSE ) {

        IntrinsicType _A, _C0, _C1, _C2, _C3;
        IntrinsicType _B[16*numElements];

        IntrinsicType _alpha;
        _alpha.fill(alpha);

        for (IndexType i=0; i<4*numElements; ++i) {
            _B[4*i  ].loadu(B+i*ldB);
            _B[4*i+1].loadu(B+i*ldB+  numElements);
            _B[4*i+2].loadu(B+i*ldB+2*numElements);
            _B[4*i+3].loadu(B+i*ldB+3*numElements);

            _B[4*i  ] = _intrinsic_mul(_alpha, _B[4*i  ]);
            _B[4*i+1] = _intrinsic_mul(_alpha, _B[4*i+1]);
            _B[4*i+2] = _intrinsic_mul(_alpha, _B[4*i+2]);
            _B[4*i+3] = _intrinsic_mul(_alpha, _B[4*i+3]);
        }

        for (IndexType i=0; i<4*numElements; ++i) {

            _C0.loadu(C+i*ldC);
            _C1.loadu(C+i*ldC+  numElements);
            _C2.loadu(C+i*ldC+2*numElements);
            _C3.loadu(C+i*ldC+3*numElements);

            for (IndexType j=0; j<4*numElements; ++j) {

                _A.fill(A[j*ldA+i]);
                _C0 = _intrinsic_add(_C0, _intrinsic_mul(_A, _B[4*j  ]));
                _C1 = _intrinsic_add(_C1, _intrinsic_mul(_A, _B[4*j+1]));
                _C2 = _intrinsic_add(_C2, _intrinsic_mul(_A, _B[4*j+2]));
                _C3 = _intrinsic_add(_C3, _intrinsic_mul(_A, _B[4*j+3]));
            }

            _C0.storeu(C+i*ldC);
            _C1.storeu(C+i*ldC+numElements);
            _C2.storeu(C+i*ldC+2*numElements);
            _C3.storeu(C+i*ldC+3*numElements);
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
gemm_real_n_n(IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const T *A, IndexType ldA,
     const T *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)

{
    CXXBLAS_DEBUG_OUT("gemv_real_n_n");

    using std::max;
    using std::min;
    using std::sqrt;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const IndexType numElements = IntrinsicType::numElements;


    const IndexType minMNK = min(m, min(n, k));
    const IndexType maxMNK = max(m, max(n, k));
    const IndexType tinyBS = 8;

    if ( minMNK==0 ) {
        return;
    }

    const IndexType N    = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL3);
    const IndexType tmpM = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL2);
    const IndexType M    = tmpM - (tmpM % numElements);

    /*
        ( C00 C01 ) = (A00 A01) ( B00 B01)
        ( C10 C11 ) = (A10 A11) ( B10 B11)
    */


    if (minMNK>=N && maxMNK>N) {

        const T one(1);

        gemm_real_n_n(N  ,  N,  N, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_n_n(N  ,  N,k-N, alpha, A+N        , ldA, B+N*ldB    , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_n_n(N  ,n-N,k-N, alpha, A+N        , ldA, B+N*(ldB+1), ldB, beta, C+N        , ldC); // C01 += A01*B11
        gemm_real_n_n(N  ,n-N,  N, alpha, A          , ldA, B+N        , ldB, one , C+N        , ldC); // C01 += A00*B01
        gemm_real_n_n(m-N,n-N,  N, alpha, A+N*ldA    , ldA, B+N        , ldB, beta, C+N*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_n_n(m-N,n-N,k-N, alpha, A+N*(ldA+1), ldA, B+N*(ldB+1), ldB, one , C+N*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_n_n(m-N,  N,k-N, alpha, A+N*(ldA+1), ldA, B+N*ldB    , ldB, beta, C+N*ldC    , ldC); // C10 += A11*B10
        gemm_real_n_n(m-N,  N,  N, alpha, A+N*ldA    , ldA, B          , ldB, one , C+N*ldC    , ldC); // C10 += A10*B00

        return;

    } else if (minMNK>=M && maxMNK>M) {

        const T one(1);

        gemm_real_n_n(M  ,  M,  M, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_n_n(M  ,  M,k-M, alpha, A+M        , ldA, B+M*ldB    , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_n_n(M  ,n-M,k-M, alpha, A+M        , ldA, B+M*(ldB+1), ldB, beta, C+M        , ldC); // C01 += A01*B11
        gemm_real_n_n(M  ,n-M,  M, alpha, A          , ldA, B+M        , ldB, one , C+M        , ldC); // C01 += A00*B01
        gemm_real_n_n(m-M,n-M,  M, alpha, A+M*ldA    , ldA, B+M        , ldB, beta, C+M*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_n_n(m-M,n-M,k-M, alpha, A+M*(ldA+1), ldA, B+M*(ldB+1), ldB, one , C+M*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_n_n(m-M,  M,k-M, alpha, A+M*(ldA+1), ldA, B+M*ldB    , ldB, beta, C+M*ldC    , ldC); // C10 += A11*B10
        gemm_real_n_n(m-M,  M,  M, alpha, A+M*ldA    , ldA, B          , ldB, one , C+M*ldC    , ldC); // C10 += A10*B00

        return;

    }


    IndexType i=0;

    for (;i+tinyBS-1<m; i+=tinyBS) {

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            scal(n, beta, C+(i+ii)*ldC, 1);
        }

        IndexType j=0;

        for (; j+tinyBS-1<n; j+=tinyBS) {
            IndexType l=0;
            for (; l+tinyBS-1<k; l+=tinyBS) {
                gemm_real_n_n_block_8x8(alpha, A+i*ldA+l, ldA, B+l*ldB+j, ldB, T(1), C+i*ldC+j, ldC);
            }

            for (; l<k; ++l) {
                for (IndexType ii=0; ii<tinyBS; ++ii) {
                    for (IndexType jj=0; jj<tinyBS; ++jj) {
                        C[(i+ii)*ldC+j+jj] += alpha*A[l+(i+ii)*ldA]*B[l*ldB+j+jj];
                    }
                }
            }
        }

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            for (IndexType l=0; l<k; ++l) {
                for (IndexType jj=j; jj<n; ++jj) {
                    C[(i+ii)*ldC+jj] += alpha*A[l+(i+ii)*ldA]*B[l*ldB+jj];
                }
            }
        }
    }

    for (; i<m; ++i) {

        scal(n, beta, C+i*ldC, 1);

        for (IndexType l=0; l<k; ++l) {

            axpy(n, alpha*A[l+i*ldA],B+l*ldB, IndexType(1) ,C+i*ldC, IndexType(1));

        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
gemm_real_t_n(IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const T *A, IndexType ldA,
     const T *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)

{
    CXXBLAS_DEBUG_OUT("gemv_real_n_n");


    using std::max;
    using std::min;
    using std::sqrt;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const IndexType numElements = IntrinsicType::numElements;


    const IndexType minMNK = min(m, min(n, k));
    const IndexType maxMNK = max(m, max(n, k));

    if ( minMNK==0 ) {
        return;
    }

    const IndexType N    = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL3);
    const IndexType tmpM = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL2);
    const IndexType M    = tmpM - (tmpM % numElements);
    const IndexType tinyBS = 8;

    /*
        ( C00 C01 ) = (A00 A01) ( B00 B01)
        ( C10 C11 ) = (A10 A11) ( B10 B11)
    */


    if (minMNK>=N && maxMNK>N) {

        const T one(1);

        gemm_real_t_n(N  ,  N,  N, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_t_n(N  ,  N,k-N, alpha, A+N*ldA    , ldA, B+N*ldB    , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_t_n(N  ,n-N,k-N, alpha, A+N*ldA    , ldA, B+N*(ldB+1), ldB, beta, C+N        , ldC); // C01 += A01*B11
        gemm_real_t_n(N  ,n-N,  N, alpha, A          , ldA, B+N        , ldB, one , C+N        , ldC); // C01 += A00*B01
        gemm_real_t_n(m-N,n-N,  N, alpha, A+N        , ldA, B+N        , ldB, beta, C+N*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_t_n(m-N,n-N,k-N, alpha, A+N*(ldA+1), ldA, B+N*(ldB+1), ldB, one , C+N*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_t_n(m-N,  N,k-N, alpha, A+N*(ldA+1), ldA, B+N*ldB    , ldB, beta, C+N*ldC    , ldC); // C10 += A11*B10
        gemm_real_t_n(m-N,  N,  N, alpha, A+N        , ldA, B          , ldB, one , C+N*ldC    , ldC); // C10 += A10*B00

        return;

    } else if (minMNK>=M && maxMNK>M) {

        const T one(1);

        gemm_real_t_n(M  ,  M,  M, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_t_n(M  ,  M,k-M, alpha, A+M*ldA    , ldA, B+M*ldB    , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_t_n(M  ,n-M,k-M, alpha, A+M*ldA    , ldA, B+M*(ldB+1), ldB, beta, C+M        , ldC); // C01 += A01*B11
        gemm_real_t_n(M  ,n-M,  M, alpha, A          , ldA, B+M        , ldB, one , C+M        , ldC); // C01 += A00*B01
        gemm_real_t_n(m-M,n-M,  M, alpha, A+M        , ldA, B+M        , ldB, beta, C+M*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_t_n(m-M,n-M,k-M, alpha, A+M*(ldA+1), ldA, B+M*(ldB+1), ldB, one , C+M*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_t_n(m-M,  M,k-M, alpha, A+M*(ldA+1), ldA, B+M*ldB    , ldB, beta, C+M*ldC    , ldC); // C10 += A11*B10
        gemm_real_t_n(m-M,  M,  M, alpha, A+M        , ldA, B          , ldB, one , C+M*ldC    , ldC); // C10 += A10*B00

        return;

    }

    IndexType i=0;

    for (;i+tinyBS-1<m; i+=tinyBS) {

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            scal(n, beta, C+(i+ii)*ldC, 1);
        }

        IndexType j=0;

        for (; j+tinyBS-1<n; j+=tinyBS) {
            IndexType l=0;
            for (; l+tinyBS-1<k; l+=tinyBS) {
                gemm_real_t_n_block_8x8(alpha, A+l*ldA+i, ldA, B+l*ldB+j, ldB, T(1), C+i*ldC+j, ldC);
            }

            for (; l<k; ++l) {
                for (IndexType ii=0; ii<tinyBS; ++ii) {
                    for (IndexType jj=0; jj<tinyBS; ++jj) {
                        C[(i+ii)*ldC+j+jj] += alpha*A[l*ldA+(i+ii)]*B[l*ldB+j+jj];
                    }
                }
            }
        }

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            for (IndexType l=0; l<k; ++l) {
                for (IndexType jj=j; jj<n; ++jj) {
                    C[(i+ii)*ldC+jj] += alpha*A[l*ldA+(i+ii)]*B[l*ldB+jj];
                }
            }
        }
    }

    for (; i<m; ++i) {

        scal(n, beta, C+i*ldC, 1);

        for (IndexType l=0; l<k; ++l) {

            axpy(n, alpha*A[i+l*ldA],B+l*ldB, IndexType(1) ,C+i*ldC, IndexType(1));

        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
gemm_real_n_t(IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const T *A, IndexType ldA,
     const T *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)

{
    CXXBLAS_DEBUG_OUT("gemv_real_n_t");

    using std::max;
    using std::min;
    using std::sqrt;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const IndexType numElements = IntrinsicType::numElements;

    const IndexType minMNK = min(m, min(n, k));
    const IndexType maxMNK = max(m, max(n, k));

    if ( minMNK==0 ) {
        return;
    }

    const IndexType N    = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL3);
    const IndexType tmpM = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL2);
    const IndexType M    = tmpM - (tmpM % numElements);

    /*
        ( C00 C01 ) = (A00 A01) ( B00 B01)
        ( C10 C11 ) = (A10 A11) ( B10 B11)
    */


    if (minMNK>=N && maxMNK>N) {

        const T one(1);

        gemm_real_n_t(N  ,  N,  N, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_n_t(N  ,  N,k-N, alpha, A+N        , ldA, B+N        , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_n_t(N  ,n-N,k-N, alpha, A+N        , ldA, B+N*(ldB+1), ldB, beta, C+N        , ldC); // C01 += A01*B11
        gemm_real_n_t(N  ,n-N,  N, alpha, A          , ldA, B+N*ldB    , ldB, one , C+N        , ldC); // C01 += A00*B01
        gemm_real_n_t(m-N,n-N,  N, alpha, A+N*ldA    , ldA, B+N*ldB    , ldB, beta, C+N*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_n_t(m-N,n-N,k-N, alpha, A+N*(ldA+1), ldA, B+N*(ldB+1), ldB, one , C+N*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_n_t(m-N,  N,k-N, alpha, A+N*(ldA+1), ldA, B+N        , ldB, beta, C+N*ldC    , ldC); // C10 += A11*B10
        gemm_real_n_t(m-N,  N,  N, alpha, A+N*ldA    , ldA, B          , ldB, one , C+N*ldC    , ldC); // C10 += A10*B00

        return;

    } else if (minMNK>=M && maxMNK>M) {

        const T one(1);

        gemm_real_n_t(M  ,  M,  M, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_n_t(M  ,  M,k-M, alpha, A+M        , ldA, B+M        , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_n_t(M  ,n-M,k-M, alpha, A+M        , ldA, B+M*(ldB+1), ldB, beta, C+M        , ldC); // C01 += A01*B11
        gemm_real_n_t(M  ,n-M,  M, alpha, A          , ldA, B+M*ldB    , ldB, one , C+M        , ldC); // C01 += A00*B01
        gemm_real_n_t(m-M,n-M,  M, alpha, A+M*ldA    , ldA, B+M*ldB    , ldB, beta, C+M*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_n_t(m-M,n-M,k-M, alpha, A+M*(ldA+1), ldA, B+M*(ldB+1), ldB, one , C+M*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_n_t(m-M,  M,k-M, alpha, A+M*(ldA+1), ldA, B+M        , ldB, beta, C+M*ldC    , ldC); // C10 += A11*B10
        gemm_real_n_t(m-M,  M,  M, alpha, A+M*ldA    , ldA, B          , ldB, one , C+M*ldC    , ldC); // C10 += A10*B00

        return;

    }

    for (IndexType i=0; i<m; ++i) {
        gemv(RowMajor, NoTrans, n, k, alpha, B, ldB, A+i*ldA, IndexType(1),
                 beta, C+i*ldC, IndexType(1));
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
gemm_real_t_t(IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const T *A, IndexType ldA,
     const T *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)

{
    CXXBLAS_DEBUG_OUT("gemv_real_t_t");

    using std::max;
    using std::min;
    using std::sqrt;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const IndexType numElements = IntrinsicType::numElements;

    const IndexType minMNK = min(m, min(n, k));
    const IndexType maxMNK = max(m, max(n, k));

    if ( minMNK==0 ) {
        return;
    }

    const IndexType N    = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL3);
    const IndexType tmpM = sqrt(BlockSize<BlasFct::GEMM, T, IndexType>::NBlockL2);
    const IndexType M    = tmpM - (tmpM % numElements);
    const IndexType tinyBS = 8;

    /*
        ( C00 C01 ) = (A00 A01) ( B00 B01)
        ( C10 C11 ) = (A10 A11) ( B10 B11)
    */


    if (minMNK>=N && maxMNK>N) {

        const T one(1);

        gemm_real_t_t(N  ,  N,  N, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_t_t(N  ,  N,k-N, alpha, A+N*ldA    , ldA, B+N        , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_t_t(N  ,n-N,k-N, alpha, A+N*ldA    , ldA, B+N*(ldB+1), ldB, beta, C+N        , ldC); // C01 += A01*B11
        gemm_real_t_t(N  ,n-N,  N, alpha, A          , ldA, B+N*ldB    , ldB, one , C+N        , ldC); // C01 += A00*B01
        gemm_real_t_t(m-N,n-N,  N, alpha, A+N        , ldA, B+N*ldB    , ldB, beta, C+N*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_t_t(m-N,n-N,k-N, alpha, A+N*(ldA+1), ldA, B+N*(ldB+1), ldB, one , C+N*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_t_t(m-N,  N,k-N, alpha, A+N*(ldA+1), ldA, B+N        , ldB, beta, C+N*ldC    , ldC); // C10 += A11*B10
        gemm_real_t_t(m-N,  N,  N, alpha, A+N        , ldA, B          , ldB, one , C+N*ldC    , ldC); // C10 += A10*B00

        return;

    } else if (minMNK>=M && maxMNK>M) {

        const T one(1);

        gemm_real_t_t(M  ,  M,  M, alpha, A          , ldA, B          , ldB, beta, C          , ldC); // C00 += A00*B00
        gemm_real_t_t(M  ,  M,k-M, alpha, A+M*ldA    , ldA, B+M        , ldB, one , C          , ldC); // C00 += A01*B10
        gemm_real_t_t(M  ,n-M,k-M, alpha, A+M*ldA    , ldA, B+M*(ldB+1), ldB, beta, C+M        , ldC); // C01 += A01*B11
        gemm_real_t_t(M  ,n-M,  M, alpha, A          , ldA, B+M*ldB    , ldB, one , C+M        , ldC); // C01 += A00*B01
        gemm_real_t_t(m-M,n-M,  M, alpha, A+M        , ldA, B+M*ldB    , ldB, beta, C+M*(ldC+1), ldC); // C11 += A10*B01
        gemm_real_t_t(m-M,n-M,k-M, alpha, A+M*(ldA+1), ldA, B+M*(ldB+1), ldB, one , C+M*(ldC+1), ldC); // C11 += A11*B11
        gemm_real_t_t(m-M,  M,k-M, alpha, A+M*(ldA+1), ldA, B+M        , ldB, beta, C+M*ldC    , ldC); // C10 += A11*B10
        gemm_real_t_t(m-M,  M,  M, alpha, A+M        , ldA, B          , ldB, one , C+M*ldC    , ldC); // C10 += A10*B00

        return;

    }

    IndexType i=0;
    T _C[8*8];

    for (;i+tinyBS-1<m; i+=tinyBS) {

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            scal(n, beta, C+(i+ii)*ldC, 1);
        }

        IndexType j=0;

        for (; j+tinyBS-1<n; j+=tinyBS) {
            IndexType l=0;

            for (IndexType p=0; p<64; ++p) {
                _C[p] = T(0);
            }
            for (; l+tinyBS-1<k; l+=tinyBS) {
                gemm_real_n_n_block_8x8(alpha,  B+l+j*ldB, ldB, A+l*ldA+i, ldA, T(1), _C, 8);
            }

            add_transpose_8x8(_C, 8, C+i*ldC+j, ldC);

            for (; l<k; ++l) {
                for (IndexType ii=0; ii<tinyBS; ++ii) {
                    for (IndexType jj=0; jj<tinyBS; ++jj) {
                        C[(i+ii)*ldC+j+jj] += alpha*A[l*ldA+(i+ii)]*B[l+(j+jj)*ldB];
                    }
                }
            }
        }

        for (IndexType ii=0; ii<tinyBS; ++ii) {
            for (IndexType l=0; l<k; ++l) {
                for (IndexType jj=j; jj<n; ++jj) {
                    C[(i+ii)*ldC+jj] += alpha*A[l*ldA+(i+ii)]*B[l+jj*ldB];
                }
            }
        }
    }

    for (; i<m; ++i) {

        scal(n, beta, C+i*ldC, 1);

        for (IndexType l=0; l<k; ++l) {

            axpy(n, alpha*A[i+l*ldA],B+l, ldB ,C+i*ldC, IndexType(1));

        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
gemm(StorageOrder order, Transpose transA, Transpose transB,
     IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const T *A, IndexType ldA,
     const T *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)
{
    CXXBLAS_DEBUG_OUT("gemm_intrinsics [real, " INTRINSIC_NAME "]");


    if ((m==0) || (n==0)) {
        return;
    }

    if (order==ColMajor) {
        gemm(RowMajor, transB, transA,
                     n, m, k, alpha,
                     B, ldB, A, ldA,
                     beta,
                     C, ldC);
        return;
    }

    if ( transA==NoTrans || transA==Conj ) {

        if ( transB==NoTrans || transB==Conj ) {

            gemm_real_n_n(m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);

        } else if ( transB==Trans || transB==ConjTrans ) {

            gemm_real_n_t(m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        }

    } else if ( transA==Trans || transA==ConjTrans ) {

        if ( transB==NoTrans || transB==Conj ) {

            gemm_real_t_n(m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);

        } else if ( transB==Trans || transB==ConjTrans ) {

            gemm_real_t_t(m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);

        }
    }

}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC
