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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/level3/gemm/packmatrix.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template<typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
kernel_gemm_2numElementsx4(IndexType k,
                           const T &alpha,
                           const T *A, IndexType ldA,
                           const T *B, IndexType ldB,
                           const T &beta,
                           T *C, IndexType ldC )
{

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    IntrinsicType  _c_00_c_10,    _c_01_c_11,    _c_02_c_12,    _c_03_c_13,
                   _c_20_c_30,    _c_21_c_31,    _c_22_c_32,    _c_23_c_33,
                   _a_0p_a_1p,    _a_2p_a_3p,
                   _b_p0,         _b_p1,         _b_p2,         _b_p3;

    _c_00_c_10.setZero();
    _c_01_c_11.setZero();
    _c_02_c_12.setZero();
    _c_03_c_13.setZero();
    _c_20_c_30.setZero();
    _c_21_c_31.setZero();
    _c_22_c_32.setZero();
    _c_23_c_33.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        _a_0p_a_1p.loadu( A );
        _a_2p_a_3p.loadu( A+numElements );
        A += 2*numElements;

        _b_p0.fill( B[0] );   /* load and duplicate */
        _b_p1.fill( B[1] );   /* load and duplicate */
        _b_p2.fill( B[2] );   /* load and duplicate */
        _b_p3.fill( B[3] );   /* load and duplicate */

        B += 4;

        /* First row and second rows */
        _c_00_c_10 = _intrinsic_add(_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _b_p0));
        _c_01_c_11 = _intrinsic_add(_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _b_p1));
        _c_02_c_12 = _intrinsic_add(_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _b_p2));
        _c_03_c_13 = _intrinsic_add(_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _b_p3));

        /* Third and fourth rows */
        _c_20_c_30 = _intrinsic_add(_c_20_c_30, _intrinsic_mul(_a_2p_a_3p, _b_p0));
        _c_21_c_31 = _intrinsic_add(_c_21_c_31, _intrinsic_mul(_a_2p_a_3p, _b_p1));
        _c_22_c_32 = _intrinsic_add(_c_22_c_32, _intrinsic_mul(_a_2p_a_3p, _b_p2));
        _c_23_c_33 = _intrinsic_add(_c_23_c_33, _intrinsic_mul(_a_2p_a_3p, _b_p3));

    }

    IntrinsicType _tmp0, _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6, _tmp7;
    IntrinsicType _alpha(alpha);

    _tmp0.loadu(C);
    _tmp1.loadu(C+ldC);
    _tmp2.loadu(C+2*ldC);
    _tmp3.loadu(C+3*ldC);

    _c_00_c_10 = _intrinsic_add(_tmp0, _intrinsic_mul(_c_00_c_10, _alpha));
    _c_01_c_11 = _intrinsic_add(_tmp1, _intrinsic_mul(_c_01_c_11, _alpha));
    _c_02_c_12 = _intrinsic_add(_tmp2, _intrinsic_mul(_c_02_c_12, _alpha));
    _c_03_c_13 = _intrinsic_add(_tmp3, _intrinsic_mul(_c_03_c_13, _alpha));


    _c_00_c_10.storeu(C);
    _c_01_c_11.storeu(C+ldC);
    _c_02_c_12.storeu(C+2*ldC);
    _c_03_c_13.storeu(C+3*ldC);


    _tmp4.loadu(C+numElements);
    _tmp5.loadu(C+numElements+ldC);
    _tmp6.loadu(C+numElements+2*ldC);
    _tmp7.loadu(C+numElements+3*ldC);

    _c_20_c_30 = _intrinsic_add(_tmp4, _intrinsic_mul(_c_20_c_30, _alpha));
    _c_21_c_31 = _intrinsic_add(_tmp5, _intrinsic_mul(_c_21_c_31, _alpha));
    _c_22_c_32 = _intrinsic_add(_tmp6, _intrinsic_mul(_c_22_c_32, _alpha));
    _c_23_c_33 = _intrinsic_add(_tmp7, _intrinsic_mul(_c_23_c_33, _alpha));

    _c_20_c_30.storeu(C+numElements);
    _c_21_c_31.storeu(C+numElements+ldC);
    _c_22_c_32.storeu(C+numElements+2*ldC);
    _c_23_c_33.storeu(C+numElements+3*ldC);

}

// Complex Kernel
template<typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
kernel_gemm_2numElementsx4(IndexType k,
                           const T &alpha,
                           const T *A, IndexType ldA,
                           const T *B, IndexType ldB,
                           const T &beta,
                           T *C, IndexType ldC )
{

    using std::real;
    using std::imag;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

    const int numElements = IntrinsicType::numElements;

    IntrinsicType  _c_00_c_10,    _c_01_c_11,    _c_02_c_12,    _c_03_c_13,
                   _c_20_c_30,    _c_21_c_31,    _c_22_c_32,    _c_23_c_33,
                   _a_0p_a_1p,    _a_2p_a_3p;
    IntrinsicPrimitiveType  _real_b_p0,   _real_b_p1,   _real_b_p2,   _real_b_p3,
                            _imag_b_p0,   _imag_b_p1,   _imag_b_p2,   _imag_b_p3;

    _c_00_c_10.setZero();
    _c_01_c_11.setZero();
    _c_02_c_12.setZero();
    _c_03_c_13.setZero();
    _c_20_c_30.setZero();
    _c_21_c_31.setZero();
    _c_22_c_32.setZero();
    _c_23_c_33.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        _a_0p_a_1p.loadu( A );
        _a_2p_a_3p.loadu( A+numElements );

        A += 2*numElements;

        _real_b_p0.fill( real(B[0]) );
        _imag_b_p0.fill( imag(B[0]) );
        _real_b_p1.fill( real(B[1]) );
        _imag_b_p1.fill( imag(B[1]) );
        _real_b_p2.fill( real(B[2]) );
        _imag_b_p2.fill( imag(B[2]) );
        _real_b_p3.fill( real(B[3]) );
        _imag_b_p3.fill( imag(B[3]) );

        B += 4;

        /* First row and second rows */
        _c_00_c_10 = _intrinsic_add(_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _real_b_p0));
        _c_01_c_11 = _intrinsic_add(_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _real_b_p1));
        _c_02_c_12 = _intrinsic_add(_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _real_b_p2));
        _c_03_c_13 = _intrinsic_add(_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _real_b_p3));

        _a_0p_a_1p = _intrinsic_swap_real_imag(_a_0p_a_1p);

        _c_00_c_10 = _intrinsic_addsub(_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _imag_b_p0));
        _c_01_c_11 = _intrinsic_addsub(_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _imag_b_p1));
        _c_02_c_12 = _intrinsic_addsub(_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _imag_b_p2));
        _c_03_c_13 = _intrinsic_addsub(_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _imag_b_p3));

        /* Third and fourth rows */
        _c_20_c_30 = _intrinsic_add(_c_20_c_30, _intrinsic_mul(_a_2p_a_3p, _real_b_p0));
        _c_21_c_31 = _intrinsic_add(_c_21_c_31, _intrinsic_mul(_a_2p_a_3p, _real_b_p1));
        _c_22_c_32 = _intrinsic_add(_c_22_c_32, _intrinsic_mul(_a_2p_a_3p, _real_b_p2));
        _c_23_c_33 = _intrinsic_add(_c_23_c_33, _intrinsic_mul(_a_2p_a_3p, _real_b_p3));

        _a_2p_a_3p = _intrinsic_swap_real_imag(_a_2p_a_3p);

        _c_20_c_30 = _intrinsic_addsub(_c_20_c_30, _intrinsic_mul(_a_2p_a_3p, _imag_b_p0));
        _c_21_c_31 = _intrinsic_addsub(_c_21_c_31, _intrinsic_mul(_a_2p_a_3p, _imag_b_p1));
        _c_22_c_32 = _intrinsic_addsub(_c_22_c_32, _intrinsic_mul(_a_2p_a_3p, _imag_b_p2));
        _c_23_c_33 = _intrinsic_addsub(_c_23_c_33, _intrinsic_mul(_a_2p_a_3p, _imag_b_p3));

    }

    IntrinsicType _tmp0, _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6, _tmp7;
    IntrinsicPrimitiveType _real_alpha(real(alpha)), _imag_alpha(imag(alpha));

    _tmp0.loadu(C);
    _tmp1.loadu(C+ldC);
    _tmp2.loadu(C+2*ldC);
    _tmp3.loadu(C+3*ldC);

    _tmp0 = _intrinsic_add(_tmp0, _intrinsic_mul(_c_00_c_10, _real_alpha));
    _tmp1 = _intrinsic_add(_tmp1, _intrinsic_mul(_c_01_c_11, _real_alpha));
    _tmp2 = _intrinsic_add(_tmp2, _intrinsic_mul(_c_02_c_12, _real_alpha));
    _tmp3 = _intrinsic_add(_tmp3, _intrinsic_mul(_c_03_c_13, _real_alpha));

    _c_00_c_10 = _intrinsic_swap_real_imag(_c_00_c_10);
    _c_01_c_11 = _intrinsic_swap_real_imag(_c_01_c_11);
    _c_02_c_12 = _intrinsic_swap_real_imag(_c_02_c_12);
    _c_03_c_13 = _intrinsic_swap_real_imag(_c_03_c_13);

    _tmp0 = _intrinsic_addsub(_tmp0, _intrinsic_mul(_c_00_c_10, _imag_alpha));
    _tmp1 = _intrinsic_addsub(_tmp1, _intrinsic_mul(_c_01_c_11, _imag_alpha));
    _tmp2 = _intrinsic_addsub(_tmp2, _intrinsic_mul(_c_02_c_12, _imag_alpha));
    _tmp3 = _intrinsic_addsub(_tmp3, _intrinsic_mul(_c_03_c_13, _imag_alpha));

    _tmp0.storeu(C);
    _tmp1.storeu(C+ldC);
    _tmp2.storeu(C+2*ldC);
    _tmp3.storeu(C+3*ldC);


    _tmp4.loadu(C+numElements);
    _tmp5.loadu(C+numElements+ldC);
    _tmp6.loadu(C+numElements+2*ldC);
    _tmp7.loadu(C+numElements+3*ldC);

    _tmp4 = _intrinsic_add(_tmp4, _intrinsic_mul(_c_20_c_30, _real_alpha));
    _tmp5 = _intrinsic_add(_tmp5, _intrinsic_mul(_c_21_c_31, _real_alpha));
    _tmp6 = _intrinsic_add(_tmp6, _intrinsic_mul(_c_22_c_32, _real_alpha));
    _tmp7 = _intrinsic_add(_tmp7, _intrinsic_mul(_c_23_c_33, _real_alpha));

    _c_20_c_30 = _intrinsic_swap_real_imag(_c_20_c_30);
    _c_21_c_31 = _intrinsic_swap_real_imag(_c_21_c_31);
    _c_22_c_32 = _intrinsic_swap_real_imag(_c_22_c_32);
    _c_23_c_33 = _intrinsic_swap_real_imag(_c_23_c_33);

    _tmp4 = _intrinsic_addsub(_tmp4, _intrinsic_mul(_c_20_c_30, _imag_alpha));
    _tmp5 = _intrinsic_addsub(_tmp5, _intrinsic_mul(_c_21_c_31, _imag_alpha));
    _tmp6 = _intrinsic_addsub(_tmp6, _intrinsic_mul(_c_22_c_32, _imag_alpha));
    _tmp7 = _intrinsic_addsub(_tmp7, _intrinsic_mul(_c_23_c_33, _imag_alpha));

    _tmp4.storeu(C+numElements);
    _tmp5.storeu(C+numElements+ldC);
    _tmp6.storeu(C+numElements+2*ldC);
    _tmp7.storeu(C+numElements+3*ldC);

}

// Mixed Kernel: A and B real, C complex
template<typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
kernel_gemm_2numElementsx4(IndexType k,
                           const std::complex<T> &alpha,
                           const T *A, IndexType ldA,
                           const T *B, IndexType ldB,
                           const std::complex<T> &beta,
                           std::complex<T> *C, IndexType ldC )
{

    using std::real;
    using std::imag;

    typedef Intrinsics<std::complex<T>, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

    const int numElements = IntrinsicType::numElements;

    IntrinsicPrimitiveType  _c_00_c_10,    _c_01_c_11,    _c_02_c_12,    _c_03_c_13,
                            _a_0p_a_1p,    _a_2p_a_3p;
    IntrinsicPrimitiveType  _real_b_p0,   _real_b_p1,   _real_b_p2,   _real_b_p3;

    _c_00_c_10.setZero();
    _c_01_c_11.setZero();
    _c_02_c_12.setZero();
    _c_03_c_13.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        _a_0p_a_1p.loadu( A );

        A += 2*numElements;

        _real_b_p0.fill( B[0] );
        _real_b_p1.fill( B[1] );
        _real_b_p2.fill( B[2] );
        _real_b_p3.fill( B[3] );

        B += 4;

        /* First row and second rows */
        _c_00_c_10 = _intrinsic_add(_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _real_b_p0));
        _c_01_c_11 = _intrinsic_add(_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _real_b_p1));
        _c_02_c_12 = _intrinsic_add(_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _real_b_p2));
        _c_03_c_13 = _intrinsic_add(_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _real_b_p3));

    }

    T *tmp = new T[4*2*numElements];

    _c_00_c_10.storeu(tmp);
    _c_01_c_11.storeu(tmp+2*numElements);
    _c_02_c_12.storeu(tmp+4*numElements);
    _c_03_c_13.storeu(tmp+6*numElements);

    for (IndexType i=0; i<4; ++i) {
        for (IndexType j=0; j<2*numElements; ++j) {
            C[i*ldC+j] += alpha*tmp[i*2*numElements+j];
        }
    }

    delete[] tmp;
}

// Mixed Kernel: B real, A and C complex
template<typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
kernel_gemm_2numElementsx4(IndexType k,
                           const std::complex<T> &alpha,
                           const std::complex<T> *A, IndexType ldA,
                           const T *B, IndexType ldB,
                           const std::complex<T> &beta,
                           std::complex<T> *C, IndexType ldC )
{

    using std::real;
    using std::imag;

    typedef Intrinsics<std::complex<T>, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

    const int numElements = IntrinsicType::numElements;

    IntrinsicType  _c_00_c_10,    _c_01_c_11,    _c_02_c_12,    _c_03_c_13,
                   _c_20_c_30,    _c_21_c_31,    _c_22_c_32,    _c_23_c_33,
                   _a_0p_a_1p,    _a_2p_a_3p;
    IntrinsicPrimitiveType  _real_b_p0,   _real_b_p1,   _real_b_p2,   _real_b_p3;

    _c_00_c_10.setZero();
    _c_01_c_11.setZero();
    _c_02_c_12.setZero();
    _c_03_c_13.setZero();
    _c_20_c_30.setZero();
    _c_21_c_31.setZero();
    _c_22_c_32.setZero();
    _c_23_c_33.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        _a_0p_a_1p.loadu( A );
        _a_2p_a_3p.loadu( A+numElements );

        A += 2*numElements;

        _real_b_p0.fill( B[0] );
        _real_b_p1.fill( B[1] );
        _real_b_p2.fill( B[2] );
        _real_b_p3.fill( B[3] );

        B += 4;

        /* First row and second rows */
        _c_00_c_10 = _intrinsic_add(_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _real_b_p0));
        _c_01_c_11 = _intrinsic_add(_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _real_b_p1));
        _c_02_c_12 = _intrinsic_add(_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _real_b_p2));
        _c_03_c_13 = _intrinsic_add(_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _real_b_p3));

        /* Third and fourth rows */
        _c_20_c_30 = _intrinsic_add(_c_20_c_30, _intrinsic_mul(_a_2p_a_3p, _real_b_p0));
        _c_21_c_31 = _intrinsic_add(_c_21_c_31, _intrinsic_mul(_a_2p_a_3p, _real_b_p1));
        _c_22_c_32 = _intrinsic_add(_c_22_c_32, _intrinsic_mul(_a_2p_a_3p, _real_b_p2));
        _c_23_c_33 = _intrinsic_add(_c_23_c_33, _intrinsic_mul(_a_2p_a_3p, _real_b_p3));

    }

    IntrinsicType _tmp0, _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6, _tmp7;
    IntrinsicPrimitiveType _real_alpha(real(alpha)), _imag_alpha(imag(alpha));

    _tmp0.loadu(C);
    _tmp1.loadu(C+ldC);
    _tmp2.loadu(C+2*ldC);
    _tmp3.loadu(C+3*ldC);

    _tmp0 = _intrinsic_add(_tmp0, _intrinsic_mul(_c_00_c_10, _real_alpha));
    _tmp1 = _intrinsic_add(_tmp1, _intrinsic_mul(_c_01_c_11, _real_alpha));
    _tmp2 = _intrinsic_add(_tmp2, _intrinsic_mul(_c_02_c_12, _real_alpha));
    _tmp3 = _intrinsic_add(_tmp3, _intrinsic_mul(_c_03_c_13, _real_alpha));

    _c_00_c_10 = _intrinsic_swap_real_imag(_c_00_c_10);
    _c_01_c_11 = _intrinsic_swap_real_imag(_c_01_c_11);
    _c_02_c_12 = _intrinsic_swap_real_imag(_c_02_c_12);
    _c_03_c_13 = _intrinsic_swap_real_imag(_c_03_c_13);

    _tmp0 = _intrinsic_addsub(_tmp0, _intrinsic_mul(_c_00_c_10, _imag_alpha));
    _tmp1 = _intrinsic_addsub(_tmp1, _intrinsic_mul(_c_01_c_11, _imag_alpha));
    _tmp2 = _intrinsic_addsub(_tmp2, _intrinsic_mul(_c_02_c_12, _imag_alpha));
    _tmp3 = _intrinsic_addsub(_tmp3, _intrinsic_mul(_c_03_c_13, _imag_alpha));

    _tmp0.storeu(C);
    _tmp1.storeu(C+ldC);
    _tmp2.storeu(C+2*ldC);
    _tmp3.storeu(C+3*ldC);


    _tmp4.loadu(C+numElements);
    _tmp5.loadu(C+numElements+ldC);
    _tmp6.loadu(C+numElements+2*ldC);
    _tmp7.loadu(C+numElements+3*ldC);

    _tmp4 = _intrinsic_add(_tmp4, _intrinsic_mul(_c_20_c_30, _real_alpha));
    _tmp5 = _intrinsic_add(_tmp5, _intrinsic_mul(_c_21_c_31, _real_alpha));
    _tmp6 = _intrinsic_add(_tmp6, _intrinsic_mul(_c_22_c_32, _real_alpha));
    _tmp7 = _intrinsic_add(_tmp7, _intrinsic_mul(_c_23_c_33, _real_alpha));

    _c_20_c_30 = _intrinsic_swap_real_imag(_c_20_c_30);
    _c_21_c_31 = _intrinsic_swap_real_imag(_c_21_c_31);
    _c_22_c_32 = _intrinsic_swap_real_imag(_c_22_c_32);
    _c_23_c_33 = _intrinsic_swap_real_imag(_c_23_c_33);

    _tmp4 = _intrinsic_addsub(_tmp4, _intrinsic_mul(_c_20_c_30, _imag_alpha));
    _tmp5 = _intrinsic_addsub(_tmp5, _intrinsic_mul(_c_21_c_31, _imag_alpha));
    _tmp6 = _intrinsic_addsub(_tmp6, _intrinsic_mul(_c_22_c_32, _imag_alpha));
    _tmp7 = _intrinsic_addsub(_tmp7, _intrinsic_mul(_c_23_c_33, _imag_alpha));

    _tmp4.storeu(C+numElements);
    _tmp5.storeu(C+numElements+ldC);
    _tmp6.storeu(C+numElements+2*ldC);
    _tmp7.storeu(C+numElements+3*ldC);

}


// Mixed Kernel: A real, B and C complex
template<typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
kernel_gemm_2numElementsx4(IndexType k,
                           const std::complex<T> &alpha,
                           const T *A, IndexType ldA,
                           const std::complex<T> *B, IndexType ldB,
                           const std::complex<T> &beta,
                           std::complex<T> *C, IndexType ldC )
{

    using std::real;
    using std::imag;

    typedef Intrinsics<std::complex<T>, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

    const int numElements = IntrinsicType::numElements;

    IntrinsicPrimitiveType  _real_c_00_c_10,    _real_c_01_c_11,    _real_c_02_c_12,    _real_c_03_c_13,
                            _imag_c_00_c_10,    _imag_c_01_c_11,    _imag_c_02_c_12,    _imag_c_03_c_13,
                             _a_0p_a_1p;
    IntrinsicPrimitiveType  _real_b_p0,   _real_b_p1,   _real_b_p2,   _real_b_p3,
                            _imag_b_p0,   _imag_b_p1,   _imag_b_p2,   _imag_b_p3;

    _real_c_00_c_10.setZero();
    _real_c_01_c_11.setZero();
    _real_c_02_c_12.setZero();
    _real_c_03_c_13.setZero();
    _imag_c_00_c_10.setZero();
    _imag_c_01_c_11.setZero();
    _imag_c_02_c_12.setZero();
    _imag_c_03_c_13.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        _a_0p_a_1p.loadu( A );
        A += 2*numElements;

        _real_b_p0.fill( real(B[0]) );
        _imag_b_p0.fill( imag(B[0]) );
        _real_b_p1.fill( real(B[1]) );
        _imag_b_p1.fill( imag(B[1]) );
        _real_b_p2.fill( real(B[2]) );
        _imag_b_p2.fill( imag(B[2]) );
        _real_b_p3.fill( real(B[3]) );
        _imag_b_p3.fill( imag(B[3]) );

        B += 4;

        /* First row and second rows */
        _real_c_00_c_10 = _intrinsic_add(_real_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _real_b_p0));
        _real_c_01_c_11 = _intrinsic_add(_real_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _real_b_p1));
        _real_c_02_c_12 = _intrinsic_add(_real_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _real_b_p2));
        _real_c_03_c_13 = _intrinsic_add(_real_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _real_b_p3));

        _imag_c_00_c_10 = _intrinsic_add(_imag_c_00_c_10, _intrinsic_mul(_a_0p_a_1p, _imag_b_p0));
        _imag_c_01_c_11 = _intrinsic_add(_imag_c_01_c_11, _intrinsic_mul(_a_0p_a_1p, _imag_b_p1));
        _imag_c_02_c_12 = _intrinsic_add(_imag_c_02_c_12, _intrinsic_mul(_a_0p_a_1p, _imag_b_p2));
        _imag_c_03_c_13 = _intrinsic_add(_imag_c_03_c_13, _intrinsic_mul(_a_0p_a_1p, _imag_b_p3));

    }
    T *real_tmp = new T[4*2*numElements];
    T *imag_tmp = new T[4*2*numElements];

    _real_c_00_c_10.storeu(real_tmp);
    _real_c_01_c_11.storeu(real_tmp+2*numElements);
    _real_c_02_c_12.storeu(real_tmp+4*numElements);
    _real_c_03_c_13.storeu(real_tmp+6*numElements);

    _imag_c_00_c_10.storeu(imag_tmp);
    _imag_c_01_c_11.storeu(imag_tmp+2*numElements);
    _imag_c_02_c_12.storeu(imag_tmp+4*numElements);
    _imag_c_03_c_13.storeu(imag_tmp+6*numElements);

    for (IndexType i=0; i<4; ++i) {
        for (IndexType j=0; j<2*numElements; ++j) {
            C[i*ldC+j] += alpha*std::complex<T>(real_tmp[i*2*numElements+j],imag_tmp[i*2*numElements+j]);
        }
    }

    delete[] real_tmp;
    delete[] imag_tmp;
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_KERNELREALGEMM_TCC
