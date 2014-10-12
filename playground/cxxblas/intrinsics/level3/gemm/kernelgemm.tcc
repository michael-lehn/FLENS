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

    IntrinsicType  c_00_c_10_,    c_01_c_11_,    c_02_c_12_,    c_03_c_13_,
                   c_20_c_30_,    c_21_c_31_,    c_22_c_32_,    c_23_c_33_,
                   a_0p_a_1p_,    a_2p_a_3p_,
                   b_p0_,         b_p1_,         b_p2_,         b_p3_;

    c_00_c_10_.setZero();
    c_01_c_11_.setZero();
    c_02_c_12_.setZero();
    c_03_c_13_.setZero();
    c_20_c_30_.setZero();
    c_21_c_31_.setZero();
    c_22_c_32_.setZero();
    c_23_c_33_.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        a_0p_a_1p_.loadu( A );
        a_2p_a_3p_.loadu( A+numElements );
        A += 2*numElements;

        b_p0_.fill( B[0] );   /* load and duplicate */
        b_p1_.fill( B[1] );   /* load and duplicate */
        b_p2_.fill( B[2] );   /* load and duplicate */
        b_p3_.fill( B[3] );   /* load and duplicate */

        B += 4;

        /* First row and second rows */
        c_00_c_10_ = intrinsic_add_(c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, b_p0_));
        c_01_c_11_ = intrinsic_add_(c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, b_p1_));
        c_02_c_12_ = intrinsic_add_(c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, b_p2_));
        c_03_c_13_ = intrinsic_add_(c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, b_p3_));

        /* Third and fourth rows */
        c_20_c_30_ = intrinsic_add_(c_20_c_30_, intrinsic_mul_(a_2p_a_3p_, b_p0_));
        c_21_c_31_ = intrinsic_add_(c_21_c_31_, intrinsic_mul_(a_2p_a_3p_, b_p1_));
        c_22_c_32_ = intrinsic_add_(c_22_c_32_, intrinsic_mul_(a_2p_a_3p_, b_p2_));
        c_23_c_33_ = intrinsic_add_(c_23_c_33_, intrinsic_mul_(a_2p_a_3p_, b_p3_));

    }

    IntrinsicType tmp0_, tmp1_, tmp2_, tmp3_, tmp4_, tmp5_, tmp6_, tmp7_;
    IntrinsicType alpha_(alpha);

    tmp0_.loadu(C);
    tmp1_.loadu(C+ldC);
    tmp2_.loadu(C+2*ldC);
    tmp3_.loadu(C+3*ldC);

    c_00_c_10_ = intrinsic_add_(tmp0_, intrinsic_mul_(c_00_c_10_, alpha_));
    c_01_c_11_ = intrinsic_add_(tmp1_, intrinsic_mul_(c_01_c_11_, alpha_));
    c_02_c_12_ = intrinsic_add_(tmp2_, intrinsic_mul_(c_02_c_12_, alpha_));
    c_03_c_13_ = intrinsic_add_(tmp3_, intrinsic_mul_(c_03_c_13_, alpha_));


    c_00_c_10_.storeu(C);
    c_01_c_11_.storeu(C+ldC);
    c_02_c_12_.storeu(C+2*ldC);
    c_03_c_13_.storeu(C+3*ldC);


    tmp4_.loadu(C+numElements);
    tmp5_.loadu(C+numElements+ldC);
    tmp6_.loadu(C+numElements+2*ldC);
    tmp7_.loadu(C+numElements+3*ldC);

    c_20_c_30_ = intrinsic_add_(tmp4_, intrinsic_mul_(c_20_c_30_, alpha_));
    c_21_c_31_ = intrinsic_add_(tmp5_, intrinsic_mul_(c_21_c_31_, alpha_));
    c_22_c_32_ = intrinsic_add_(tmp6_, intrinsic_mul_(c_22_c_32_, alpha_));
    c_23_c_33_ = intrinsic_add_(tmp7_, intrinsic_mul_(c_23_c_33_, alpha_));

    c_20_c_30_.storeu(C+numElements);
    c_21_c_31_.storeu(C+numElements+ldC);
    c_22_c_32_.storeu(C+numElements+2*ldC);
    c_23_c_33_.storeu(C+numElements+3*ldC);

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

    IntrinsicType  c_00_c_10_,    c_01_c_11_,    c_02_c_12_,    c_03_c_13_,
                   c_20_c_30_,    c_21_c_31_,    c_22_c_32_,    c_23_c_33_,
                   a_0p_a_1p_,    a_2p_a_3p_;
    IntrinsicPrimitiveType  real_b_p0_,   real_b_p1_,   real_b_p2_,   real_b_p3_,
                            imag_b_p0_,   imag_b_p1_,   imag_b_p2_,   imag_b_p3_;

    c_00_c_10_.setZero();
    c_01_c_11_.setZero();
    c_02_c_12_.setZero();
    c_03_c_13_.setZero();
    c_20_c_30_.setZero();
    c_21_c_31_.setZero();
    c_22_c_32_.setZero();
    c_23_c_33_.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        a_0p_a_1p_.loadu( A );
        a_2p_a_3p_.loadu( A+numElements );

        A += 2*numElements;

        real_b_p0_.fill( real(B[0]) );
        imag_b_p0_.fill( imag(B[0]) );
        real_b_p1_.fill( real(B[1]) );
        imag_b_p1_.fill( imag(B[1]) );
        real_b_p2_.fill( real(B[2]) );
        imag_b_p2_.fill( imag(B[2]) );
        real_b_p3_.fill( real(B[3]) );
        imag_b_p3_.fill( imag(B[3]) );

        B += 4;

        /* First row and second rows */
        c_00_c_10_ = intrinsic_add_(c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, real_b_p0_));
        c_01_c_11_ = intrinsic_add_(c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, real_b_p1_));
        c_02_c_12_ = intrinsic_add_(c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, real_b_p2_));
        c_03_c_13_ = intrinsic_add_(c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, real_b_p3_));

        a_0p_a_1p_ = intrinsic_swap_real_imag_(a_0p_a_1p_);

        c_00_c_10_ = intrinsic_addsub_(c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, imag_b_p0_));
        c_01_c_11_ = intrinsic_addsub_(c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, imag_b_p1_));
        c_02_c_12_ = intrinsic_addsub_(c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, imag_b_p2_));
        c_03_c_13_ = intrinsic_addsub_(c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, imag_b_p3_));

        /* Third and fourth rows */
        c_20_c_30_ = intrinsic_add_(c_20_c_30_, intrinsic_mul_(a_2p_a_3p_, real_b_p0_));
        c_21_c_31_ = intrinsic_add_(c_21_c_31_, intrinsic_mul_(a_2p_a_3p_, real_b_p1_));
        c_22_c_32_ = intrinsic_add_(c_22_c_32_, intrinsic_mul_(a_2p_a_3p_, real_b_p2_));
        c_23_c_33_ = intrinsic_add_(c_23_c_33_, intrinsic_mul_(a_2p_a_3p_, real_b_p3_));

        a_2p_a_3p_ = intrinsic_swap_real_imag_(a_2p_a_3p_);

        c_20_c_30_ = intrinsic_addsub_(c_20_c_30_, intrinsic_mul_(a_2p_a_3p_, imag_b_p0_));
        c_21_c_31_ = intrinsic_addsub_(c_21_c_31_, intrinsic_mul_(a_2p_a_3p_, imag_b_p1_));
        c_22_c_32_ = intrinsic_addsub_(c_22_c_32_, intrinsic_mul_(a_2p_a_3p_, imag_b_p2_));
        c_23_c_33_ = intrinsic_addsub_(c_23_c_33_, intrinsic_mul_(a_2p_a_3p_, imag_b_p3_));

    }

    IntrinsicType tmp0_, tmp1_, tmp2_, tmp3_, tmp4_, tmp5_, tmp6_, tmp7_;
    IntrinsicPrimitiveType real_alpha_(real(alpha)), imag_alpha_(imag(alpha));

    tmp0_.loadu(C);
    tmp1_.loadu(C+ldC);
    tmp2_.loadu(C+2*ldC);
    tmp3_.loadu(C+3*ldC);

    tmp0_ = intrinsic_add_(tmp0_, intrinsic_mul_(c_00_c_10_, real_alpha_));
    tmp1_ = intrinsic_add_(tmp1_, intrinsic_mul_(c_01_c_11_, real_alpha_));
    tmp2_ = intrinsic_add_(tmp2_, intrinsic_mul_(c_02_c_12_, real_alpha_));
    tmp3_ = intrinsic_add_(tmp3_, intrinsic_mul_(c_03_c_13_, real_alpha_));

    c_00_c_10_ = intrinsic_swap_real_imag_(c_00_c_10_);
    c_01_c_11_ = intrinsic_swap_real_imag_(c_01_c_11_);
    c_02_c_12_ = intrinsic_swap_real_imag_(c_02_c_12_);
    c_03_c_13_ = intrinsic_swap_real_imag_(c_03_c_13_);

    tmp0_ = intrinsic_addsub_(tmp0_, intrinsic_mul_(c_00_c_10_, imag_alpha_));
    tmp1_ = intrinsic_addsub_(tmp1_, intrinsic_mul_(c_01_c_11_, imag_alpha_));
    tmp2_ = intrinsic_addsub_(tmp2_, intrinsic_mul_(c_02_c_12_, imag_alpha_));
    tmp3_ = intrinsic_addsub_(tmp3_, intrinsic_mul_(c_03_c_13_, imag_alpha_));

    tmp0_.storeu(C);
    tmp1_.storeu(C+ldC);
    tmp2_.storeu(C+2*ldC);
    tmp3_.storeu(C+3*ldC);


    tmp4_.loadu(C+numElements);
    tmp5_.loadu(C+numElements+ldC);
    tmp6_.loadu(C+numElements+2*ldC);
    tmp7_.loadu(C+numElements+3*ldC);

    tmp4_ = intrinsic_add_(tmp4_, intrinsic_mul_(c_20_c_30_, real_alpha_));
    tmp5_ = intrinsic_add_(tmp5_, intrinsic_mul_(c_21_c_31_, real_alpha_));
    tmp6_ = intrinsic_add_(tmp6_, intrinsic_mul_(c_22_c_32_, real_alpha_));
    tmp7_ = intrinsic_add_(tmp7_, intrinsic_mul_(c_23_c_33_, real_alpha_));

    c_20_c_30_ = intrinsic_swap_real_imag_(c_20_c_30_);
    c_21_c_31_ = intrinsic_swap_real_imag_(c_21_c_31_);
    c_22_c_32_ = intrinsic_swap_real_imag_(c_22_c_32_);
    c_23_c_33_ = intrinsic_swap_real_imag_(c_23_c_33_);

    tmp4_ = intrinsic_addsub_(tmp4_, intrinsic_mul_(c_20_c_30_, imag_alpha_));
    tmp5_ = intrinsic_addsub_(tmp5_, intrinsic_mul_(c_21_c_31_, imag_alpha_));
    tmp6_ = intrinsic_addsub_(tmp6_, intrinsic_mul_(c_22_c_32_, imag_alpha_));
    tmp7_ = intrinsic_addsub_(tmp7_, intrinsic_mul_(c_23_c_33_, imag_alpha_));

    tmp4_.storeu(C+numElements);
    tmp5_.storeu(C+numElements+ldC);
    tmp6_.storeu(C+numElements+2*ldC);
    tmp7_.storeu(C+numElements+3*ldC);

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

    IntrinsicPrimitiveType  c_00_c_10_,    c_01_c_11_,    c_02_c_12_,    c_03_c_13_,
                            a_0p_a_1p_,    a_2p_a_3p_;
    IntrinsicPrimitiveType  real_b_p0_,   real_b_p1_,   real_b_p2_,   real_b_p3_;

    c_00_c_10_.setZero();
    c_01_c_11_.setZero();
    c_02_c_12_.setZero();
    c_03_c_13_.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        a_0p_a_1p_.loadu( A );

        A += 2*numElements;

        real_b_p0_.fill( B[0] );
        real_b_p1_.fill( B[1] );
        real_b_p2_.fill( B[2] );
        real_b_p3_.fill( B[3] );

        B += 4;

        /* First row and second rows */
        c_00_c_10_ = intrinsic_add_(c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, real_b_p0_));
        c_01_c_11_ = intrinsic_add_(c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, real_b_p1_));
        c_02_c_12_ = intrinsic_add_(c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, real_b_p2_));
        c_03_c_13_ = intrinsic_add_(c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, real_b_p3_));

    }

    T *tmp = new T[4*2*numElements];

    c_00_c_10_.storeu(tmp);
    c_01_c_11_.storeu(tmp+2*numElements);
    c_02_c_12_.storeu(tmp+4*numElements);
    c_03_c_13_.storeu(tmp+6*numElements);

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

    IntrinsicType  c_00_c_10_,    c_01_c_11_,    c_02_c_12_,    c_03_c_13_,
                   c_20_c_30_,    c_21_c_31_,    c_22_c_32_,    c_23_c_33_,
                   a_0p_a_1p_,    a_2p_a_3p_;
    IntrinsicPrimitiveType  real_b_p0_,   real_b_p1_,   real_b_p2_,   real_b_p3_;

    c_00_c_10_.setZero();
    c_01_c_11_.setZero();
    c_02_c_12_.setZero();
    c_03_c_13_.setZero();
    c_20_c_30_.setZero();
    c_21_c_31_.setZero();
    c_22_c_32_.setZero();
    c_23_c_33_.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        a_0p_a_1p_.loadu( A );
        a_2p_a_3p_.loadu( A+numElements );

        A += 2*numElements;

        real_b_p0_.fill( B[0] );
        real_b_p1_.fill( B[1] );
        real_b_p2_.fill( B[2] );
        real_b_p3_.fill( B[3] );

        B += 4;

        /* First row and second rows */
        c_00_c_10_ = intrinsic_add_(c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, real_b_p0_));
        c_01_c_11_ = intrinsic_add_(c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, real_b_p1_));
        c_02_c_12_ = intrinsic_add_(c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, real_b_p2_));
        c_03_c_13_ = intrinsic_add_(c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, real_b_p3_));

        /* Third and fourth rows */
        c_20_c_30_ = intrinsic_add_(c_20_c_30_, intrinsic_mul_(a_2p_a_3p_, real_b_p0_));
        c_21_c_31_ = intrinsic_add_(c_21_c_31_, intrinsic_mul_(a_2p_a_3p_, real_b_p1_));
        c_22_c_32_ = intrinsic_add_(c_22_c_32_, intrinsic_mul_(a_2p_a_3p_, real_b_p2_));
        c_23_c_33_ = intrinsic_add_(c_23_c_33_, intrinsic_mul_(a_2p_a_3p_, real_b_p3_));

    }

    IntrinsicType tmp0_, tmp1_, tmp2_, tmp3_, tmp4_, tmp5_, tmp6_, tmp7_;
    IntrinsicPrimitiveType real_alpha_(real(alpha)), imag_alpha_(imag(alpha));

    tmp0_.loadu(C);
    tmp1_.loadu(C+ldC);
    tmp2_.loadu(C+2*ldC);
    tmp3_.loadu(C+3*ldC);

    tmp0_ = intrinsic_add_(tmp0_, intrinsic_mul_(c_00_c_10_, real_alpha_));
    tmp1_ = intrinsic_add_(tmp1_, intrinsic_mul_(c_01_c_11_, real_alpha_));
    tmp2_ = intrinsic_add_(tmp2_, intrinsic_mul_(c_02_c_12_, real_alpha_));
    tmp3_ = intrinsic_add_(tmp3_, intrinsic_mul_(c_03_c_13_, real_alpha_));

    c_00_c_10_ = intrinsic_swap_real_imag_(c_00_c_10_);
    c_01_c_11_ = intrinsic_swap_real_imag_(c_01_c_11_);
    c_02_c_12_ = intrinsic_swap_real_imag_(c_02_c_12_);
    c_03_c_13_ = intrinsic_swap_real_imag_(c_03_c_13_);

    tmp0_ = intrinsic_addsub_(tmp0_, intrinsic_mul_(c_00_c_10_, imag_alpha_));
    tmp1_ = intrinsic_addsub_(tmp1_, intrinsic_mul_(c_01_c_11_, imag_alpha_));
    tmp2_ = intrinsic_addsub_(tmp2_, intrinsic_mul_(c_02_c_12_, imag_alpha_));
    tmp3_ = intrinsic_addsub_(tmp3_, intrinsic_mul_(c_03_c_13_, imag_alpha_));

    tmp0_.storeu(C);
    tmp1_.storeu(C+ldC);
    tmp2_.storeu(C+2*ldC);
    tmp3_.storeu(C+3*ldC);


    tmp4_.loadu(C+numElements);
    tmp5_.loadu(C+numElements+ldC);
    tmp6_.loadu(C+numElements+2*ldC);
    tmp7_.loadu(C+numElements+3*ldC);

    tmp4_ = intrinsic_add_(tmp4_, intrinsic_mul_(c_20_c_30_, real_alpha_));
    tmp5_ = intrinsic_add_(tmp5_, intrinsic_mul_(c_21_c_31_, real_alpha_));
    tmp6_ = intrinsic_add_(tmp6_, intrinsic_mul_(c_22_c_32_, real_alpha_));
    tmp7_ = intrinsic_add_(tmp7_, intrinsic_mul_(c_23_c_33_, real_alpha_));

    c_20_c_30_ = intrinsic_swap_real_imag_(c_20_c_30_);
    c_21_c_31_ = intrinsic_swap_real_imag_(c_21_c_31_);
    c_22_c_32_ = intrinsic_swap_real_imag_(c_22_c_32_);
    c_23_c_33_ = intrinsic_swap_real_imag_(c_23_c_33_);

    tmp4_ = intrinsic_addsub_(tmp4_, intrinsic_mul_(c_20_c_30_, imag_alpha_));
    tmp5_ = intrinsic_addsub_(tmp5_, intrinsic_mul_(c_21_c_31_, imag_alpha_));
    tmp6_ = intrinsic_addsub_(tmp6_, intrinsic_mul_(c_22_c_32_, imag_alpha_));
    tmp7_ = intrinsic_addsub_(tmp7_, intrinsic_mul_(c_23_c_33_, imag_alpha_));

    tmp4_.storeu(C+numElements);
    tmp5_.storeu(C+numElements+ldC);
    tmp6_.storeu(C+numElements+2*ldC);
    tmp7_.storeu(C+numElements+3*ldC);

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

    IntrinsicPrimitiveType  real_c_00_c_10_,    real_c_01_c_11_,    real_c_02_c_12_,    real_c_03_c_13_,
                            imag_c_00_c_10_,    imag_c_01_c_11_,    imag_c_02_c_12_,    imag_c_03_c_13_,
                             a_0p_a_1p_;
    IntrinsicPrimitiveType  real_b_p0_,   real_b_p1_,   real_b_p2_,   real_b_p3_,
                            imag_b_p0_,   imag_b_p1_,   imag_b_p2_,   imag_b_p3_;

    real_c_00_c_10_.setZero();
    real_c_01_c_11_.setZero();
    real_c_02_c_12_.setZero();
    real_c_03_c_13_.setZero();
    imag_c_00_c_10_.setZero();
    imag_c_01_c_11_.setZero();
    imag_c_02_c_12_.setZero();
    imag_c_03_c_13_.setZero();

    IndexType p=0;
    for ( ; p<k; p++ ){
        a_0p_a_1p_.loadu( A );
        A += 2*numElements;

        real_b_p0_.fill( real(B[0]) );
        imag_b_p0_.fill( imag(B[0]) );
        real_b_p1_.fill( real(B[1]) );
        imag_b_p1_.fill( imag(B[1]) );
        real_b_p2_.fill( real(B[2]) );
        imag_b_p2_.fill( imag(B[2]) );
        real_b_p3_.fill( real(B[3]) );
        imag_b_p3_.fill( imag(B[3]) );

        B += 4;

        /* First row and second rows */
        real_c_00_c_10_ = intrinsic_add_(real_c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, real_b_p0_));
        real_c_01_c_11_ = intrinsic_add_(real_c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, real_b_p1_));
        real_c_02_c_12_ = intrinsic_add_(real_c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, real_b_p2_));
        real_c_03_c_13_ = intrinsic_add_(real_c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, real_b_p3_));

        imag_c_00_c_10_ = intrinsic_add_(imag_c_00_c_10_, intrinsic_mul_(a_0p_a_1p_, imag_b_p0_));
        imag_c_01_c_11_ = intrinsic_add_(imag_c_01_c_11_, intrinsic_mul_(a_0p_a_1p_, imag_b_p1_));
        imag_c_02_c_12_ = intrinsic_add_(imag_c_02_c_12_, intrinsic_mul_(a_0p_a_1p_, imag_b_p2_));
        imag_c_03_c_13_ = intrinsic_add_(imag_c_03_c_13_, intrinsic_mul_(a_0p_a_1p_, imag_b_p3_));

    }
    T *real_tmp = new T[4*2*numElements];
    T *imag_tmp = new T[4*2*numElements];

    real_c_00_c_10_.storeu(real_tmp);
    real_c_01_c_11_.storeu(real_tmp+2*numElements);
    real_c_02_c_12_.storeu(real_tmp+4*numElements);
    real_c_03_c_13_.storeu(real_tmp+6*numElements);

    imag_c_00_c_10_.storeu(imag_tmp);
    imag_c_01_c_11_.storeu(imag_tmp+2*numElements);
    imag_c_02_c_12_.storeu(imag_tmp+4*numElements);
    imag_c_03_c_13_.storeu(imag_tmp+6*numElements);

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
