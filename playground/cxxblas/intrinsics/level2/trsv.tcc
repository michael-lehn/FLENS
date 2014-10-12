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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRSV_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRSV_TCC 1

#include <cxxstd/iostream.h>

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trsv_real_up_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_up_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {
        trsv_real_up_n(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_real_n(n-NBL3, NBL3, T(-1), A+(n-NBL3), ldA, x+n-NBL3, 1, T(1),
                    x, 1);
        trsv_real_up_n(diag, n-NBL3, A, ldA, x, 1);
        return;
    } else if (n>NBL2 ) {
        trsv_real_up_n(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_real_n(n-NBL2, NBL2, T(-1), A+(n-NBL2), ldA, x+n-NBL2, 1, T(1),
                    x, 1);
        trsv_real_up_n(diag, n-NBL2, A, ldA, x, 1);
        return;
    } else if (n>NBL1 ) {
        trsv_real_up_n(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_real_n(n-NBL1, NBL1, T(-1), A+(n-NBL1), ldA, x+n-NBL1, 1, T(1),
                    x, 1);
        trsv_real_up_n(diag, n-NBL1, A, ldA, x, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=n-1; i>=0; --i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType result_;
        result_.setZero();

        T result[numElements];

        for (; j+4*numElements-1<n-i-1; j+=4*numElements) {
            A0_.loadu(A+i*(ldA+1)+1+j);
            A1_.loadu(A+i*(ldA+1)+1+j+numElements);
            A2_.loadu(A+i*(ldA+1)+1+j+2*numElements);
            A3_.loadu(A+i*(ldA+1)+1+j+3*numElements);
            x0_.loadu(x+i+j+1);
            x1_.loadu(x+i+j+1+numElements);
            x2_.loadu(x+i+j+1+2*numElements);
            x3_.loadu(x+i+j+1+3*numElements);
            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, x3_));
        }

        for (; j+numElements-1<n-i-1; j+=numElements) {
            A0_.loadu(A+i*(ldA+1)+1+j);
            x0_.loadu(x+i+j+1);
            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }

        for (; j<n-i-1; ++j) {
            x[i] -= A[i*(ldA+1)+1+j]*x[i+1+j];
        }
        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trsv_real_up_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_up_t");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_real_up_t(diag, NBL3, A, ldA, x, 1);
        gemv_real_t(NBL3, n-NBL3, T(-1), A+NBL3, ldA, x, 1, T(1), x+NBL3, 1);
        trsv_real_up_t(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;

    } else if (n>NBL2 ) {

        trsv_real_up_t(diag, NBL2, A, ldA, x, 1);
        gemv_real_t(NBL2, n-NBL2, T(-1), A+NBL2, ldA, x, 1, T(1), x+NBL2, 1);
        trsv_real_up_t(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;

    } else if (n>NBL1 ) {

        trsv_real_up_t(diag, NBL1, A, ldA, x, 1);
        gemv_real_t(NBL1, n-NBL1, T(-1), A+NBL1, ldA, x, 1, T(1), x+NBL1, 1);
        trsv_real_up_t(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;

        x_.loadu(x+i);

        IndexType j=0;

        for (; j+3<i; j+=4) {

            x0_.fill(x[j]);
            x1_.fill(x[j+1]);
            x2_.fill(x[j+2]);
            x3_.fill(x[j+3]);

            A0_.loadu(A+i+j*ldA);
            A1_.loadu(A+i+(j+1)*ldA);
            A2_.loadu(A+i+(j+2)*ldA);
            A3_.loadu(A+i+(j+3)*ldA);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, x3_));
        }

        for (; j<i; ++j) {

            x0_.fill(x[j]);

            A0_.loadu(A+i+j*ldA);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, x0_));
        }

        x_.storeu(x+i);


        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=i; j<i+ii; ++j) {
                x[i+ii] -= A[(i+ii)+j*ldA]*x[j];
            }

            if (diag==NonUnit) {
                x[i+ii] /= A[(i+ii)*(ldA+1)];
            }
        }
    }

    for (; i<n; ++i) {

        for (IndexType j=0; j<i; ++j) {
            x[i] -= A[i+j*ldA]*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trsv_real_lo_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_lo_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_real_lo_n(diag, NBL3, A, ldA, x, 1);
        gemv_real_n(n-NBL3, NBL3, T(-1), A+NBL3*ldA, ldA, x, 1, T(1),
                    x+NBL3, 1);
        trsv_real_lo_n(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;

    } else if (n>NBL2 ) {

        trsv_real_lo_n(diag, NBL2, A, ldA, x, 1);
        gemv_real_n(n-NBL2, NBL2, T(-1), A+NBL2*ldA, ldA, x, 1, T(1),
                    x+NBL2, 1);
        trsv_real_lo_n(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;

    } else if (n>NBL1 ) {

        trsv_real_lo_n(diag, NBL1, A, ldA, x, 1);
        gemv_real_n(n-NBL1, NBL1, T(-1), A+NBL1*ldA, ldA, x, 1, T(1),
                    x+NBL1, 1);
        trsv_real_lo_n(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=0; i<n; ++i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType result_;
        result_.setZero();

        T result[numElements];

        for (; j+4*numElements-1<i; j+=4*numElements) {

            A0_.loadu(A+i*ldA+j);
            A1_.loadu(A+i*ldA+j+numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

            x0_.loadu(x+j);
            x1_.loadu(x+j+numElements);
            x2_.loadu(x+j+2*numElements);
            x3_.loadu(x+j+3*numElements);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, x3_));
        }

        for (; j+numElements-1<i; j+=numElements) {
            A0_.loadu(A+i*ldA+j);

            x0_.loadu(x+j);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }
        for (; j<i; ++j) {
            x[i] -= A[i*ldA+j]*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trsv_real_lo_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_lo_t");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();


    if (n>NBL3 ) {
        trsv_real_lo_t(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_real_t(NBL3, n-NBL3, T(-1), A+(n-NBL3)*ldA, ldA, x+n-NBL3, 1, T(1),
                    x, 1);
        trsv_real_lo_t(diag, n-NBL3, A, ldA, x, 1);
        return;
    } else if (n>NBL2 ) {
        trsv_real_lo_t(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_real_t(NBL2, n-NBL2, T(-1), A+(n-NBL2)*ldA, ldA, x+n-NBL2, 1, T(1),
                    x, 1);
        trsv_real_lo_t(diag, n-NBL2, A, ldA, x, 1);
        return;
    } else if (n>NBL1 ) {
        trsv_real_lo_t(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_real_t(NBL1, n-NBL1, T(-1), A+(n-NBL1)*ldA, ldA, x+n-NBL1, 1, T(1),
                    x, 1);
        trsv_real_lo_t(diag, n-NBL1, A, ldA, x, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=n-1;

    for (; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;

        x_.loadu(x+i-numElements+1);

        IndexType j=0;

        for (; j+3<n-i-1; j+=4) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);
            A1_.loadu(A+(i+j+2)*ldA+i-numElements+1);
            A2_.loadu(A+(i+j+3)*ldA+i-numElements+1);
            A3_.loadu(A+(i+j+4)*ldA+i-numElements+1);

            x0_.fill(x[i+j+1]);
            x1_.fill(x[i+j+2]);
            x2_.fill(x[i+j+3]);
            x3_.fill(x[i+j+4]);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, x3_));
        }

        for (; j<n-i-1; ++j) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);

            x0_.fill(x[i+j+1]);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, x0_));
        }
        x_.storeu(x+i-numElements+1);

        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=0; j<ii; ++j) {
                x[i-ii] -= A[(i-ii)*(ldA+1)+(j+1)*ldA]*x[(i-ii)+j+1];
            }

            if (diag==NonUnit) {
                x[i-ii] /= A[(i-ii)*(ldA+1)];
            }
        }
    }

    for (; i>=0; --i) {

        for (IndexType j=0; j<n-i-1; ++j) {
            x[i] -= A[i*(ldA+1)+(j+1)*ldA]*x[i+j+1];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_up_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_up_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {
        trsv_complex_up_n(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_complex_n(n-NBL3, NBL3, T(-1), A+(n-NBL3), ldA, x+n-NBL3, 1, T(1),
                       x, 1);
        trsv_complex_up_n(diag, n-NBL3, A, ldA, x, 1);
        return;
    } else if (n>NBL2 ) {
        trsv_complex_up_n(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_complex_n(n-NBL2, NBL2, T(-1), A+(n-NBL2), ldA, x+n-NBL2, 1, T(1),
                       x, 1);
        trsv_complex_up_n(diag, n-NBL2, A, ldA, x, 1);
        return;
    } else if (n>NBL1 ) {
        trsv_complex_up_n(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_complex_n(n-NBL1, NBL1, T(-1), A+(n-NBL1), ldA, x+n-NBL1, 1, T(1),
                       x, 1);
        trsv_complex_up_n(diag, n-NBL1, A, ldA, x, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=n-1; i>=0; --i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;


        IntrinsicType result_;
        result_.setZero();

        T result[numElements];

        for (; j+4*numElements-1<n-i-1; j+=4*numElements) {

            A0_.loadu(A+i*(ldA+1)+1+j);
            A1_.loadu(A+i*(ldA+1)+1+j+numElements);
            A2_.loadu(A+i*(ldA+1)+1+j+2*numElements);
            A3_.loadu(A+i*(ldA+1)+1+j+3*numElements);

            x0_.loadu(x+i+j+1);
            x1_.loadu(x+i+j+1+numElements);
            x2_.loadu(x+i+j+1+2*numElements);
            x3_.loadu(x+i+j+1+3*numElements);

            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);
            real_x1_ = intrinsic_real_(x1_);
            imag_x1_ = intrinsic_imag_(x1_);
            real_x2_ = intrinsic_real_(x2_);
            imag_x2_ = intrinsic_imag_(x2_);
            real_x3_ = intrinsic_real_(x3_);
            imag_x3_ = intrinsic_imag_(x3_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, real_x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, real_x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A1_, imag_x1_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A2_, imag_x2_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A3_, imag_x3_));

        }

        for (; j+numElements-1<n-i-1; j+=numElements) {

            A0_.loadu(A+i*(ldA+1)+1+j);

            x0_.loadu(x+i+j+1);

            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }

        for (; j<n-i-1; ++j) {
            x[i] -= A[i*(ldA+1)+1+j]*x[i+1+j];
        }
        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_up_c(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_complex_c");

    using std::conj;

    ASSERT( incX==1 );
    CXXBLAS_DEBUG_OUT("trsv_complex_up_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {
        trsv_complex_up_c(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_complex_c(n-NBL3, NBL3, T(-1), A+(n-NBL3), ldA, x+n-NBL3, 1, T(1),
                       x, 1);
        trsv_complex_up_c(diag, n-NBL3, A, ldA, x, 1);
        return;
    } else if (n>NBL2 ) {
        trsv_complex_up_c(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_complex_c(n-NBL2, NBL2, T(-1), A+(n-NBL2), ldA, x+n-NBL2, 1, T(1),
                       x, 1);
        trsv_complex_up_c(diag, n-NBL2, A, ldA, x, 1);
        return;
    } else if (n>NBL1 ) {
        trsv_complex_up_c(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_complex_c(n-NBL1, NBL1, T(-1), A+(n-NBL1), ldA, x+n-NBL1, 1, T(1),
                       x, 1);
        trsv_complex_up_c(diag, n-NBL1, A, ldA, x, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=n-1; i>=0; --i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;


        IntrinsicType result_;
        result_.setZero();

        T result[numElements];

        IntrinsicPrimitiveType minusOne_;
        PT minusOne(-1);
        minusOne_.fill(minusOne);

        for (; j+4*numElements-1<n-i-1; j+=4*numElements) {

            A0_.loadu(A+i*(ldA+1)+1+j);
            A1_.loadu(A+i*(ldA+1)+1+j+numElements);
            A2_.loadu(A+i*(ldA+1)+1+j+2*numElements);
            A3_.loadu(A+i*(ldA+1)+1+j+3*numElements);

            x0_.loadu(x+i+j+1);
            x1_.loadu(x+i+j+1+numElements);
            x2_.loadu(x+i+j+1+2*numElements);
            x3_.loadu(x+i+j+1+3*numElements);

            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);
            real_x1_ = intrinsic_mul_(minusOne_, intrinsic_real_(x1_));
            imag_x1_ = intrinsic_imag_(x1_);
            real_x2_ = intrinsic_mul_(minusOne_, intrinsic_real_(x2_));
            imag_x2_ = intrinsic_imag_(x2_);
            real_x3_ = intrinsic_mul_(minusOne_, intrinsic_real_(x3_));
            imag_x3_ = intrinsic_imag_(x3_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A1_, real_x1_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A2_, real_x2_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, imag_x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, imag_x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, imag_x3_));

        }

        for (; j+numElements-1<n-i-1; j+=numElements) {

            A0_.loadu(A+i*(ldA+1)+1+j);

            x0_.loadu(x+i+j+1);

            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }

        for (; j<n-i-1; ++j) {
            x[i] -= conj(A[i*(ldA+1)+1+j])*x[i+1+j];
        }
        if (diag==NonUnit) {
            x[i] /= conj(A[i*(ldA+1)]);
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_up_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_up_t");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_complex_up_t(diag, NBL3, A, ldA, x, 1);
        gemv_complex_t(NBL3, n-NBL3, T(-1), A+NBL3, ldA, x, 1, T(1), x+NBL3, 1);
        trsv_complex_up_t(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;
    } else if (n>NBL2 ) {

        trsv_complex_up_t(diag, NBL2, A, ldA, x, 1);
        gemv_complex_t(NBL2, n-NBL2, T(-1), A+NBL2, ldA, x, 1, T(1), x+NBL2, 1);
        trsv_complex_up_t(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;
    } else if (n>NBL1 ) {

        trsv_complex_up_t(diag, NBL1, A, ldA, x, 1);
        gemv_complex_t(NBL1, n-NBL1, T(-1), A+NBL1, ldA, x, 1, T(1), x+NBL1, 1);
        trsv_complex_up_t(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i);

        IndexType j=0;

        for (; j+3<i; j+=4) {

            real_x0_.fill( real(x[j]));
            imag_x0_.fill(-imag(x[j]));
            real_x1_.fill( real(x[j+1]));
            imag_x1_.fill(-imag(x[j+1]));
            real_x2_.fill( real(x[j+2]));
            imag_x2_.fill(-imag(x[j+2]));
            real_x3_.fill( real(x[j+3]));
            imag_x3_.fill(-imag(x[j+3]));

            A0_.loadu(A+i+j*ldA);
            A1_.loadu(A+i+(j+1)*ldA);
            A2_.loadu(A+i+(j+2)*ldA);
            A3_.loadu(A+i+(j+3)*ldA);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j<i; ++j) {

            real_x0_.fill( real(x[j]));
            imag_x0_.fill(-imag(x[j]));

            A0_.loadu(A+i+j*ldA);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
        }

        x_.storeu(x+i);


        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=i; j<i+ii; ++j) {
                x[i+ii] -= A[(i+ii)+j*ldA]*x[j];
            }

            if (diag==NonUnit) {
                x[i+ii] /= A[(i+ii)*(ldA+1)];
            }
        }
    }

    for (; i<n; ++i) {

        for (IndexType j=0; j<i; ++j) {
            x[i] -= A[i+j*ldA]*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_up_ct(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_up_ct");

    using std::conj;
    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_complex_up_ct(diag, NBL3, A, ldA, x, 1);
        gemv_complex_ct(NBL3, n-NBL3, T(-1), A+NBL3, ldA, x, 1, T(1),
                        x+NBL3, 1);
        trsv_complex_up_ct(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;
    } else  if (n>NBL2 ) {

        trsv_complex_up_ct(diag, NBL2, A, ldA, x, 1);
        gemv_complex_ct(NBL2, n-NBL2, T(-1), A+NBL2, ldA, x, 1, T(1),
                        x+NBL2, 1);
        trsv_complex_up_ct(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;
    } else if (n>NBL1 ) {

        trsv_complex_up_ct(diag, NBL1, A, ldA, x, 1);
        gemv_complex_ct(NBL1, n-NBL1, T(-1), A+NBL1, ldA, x, 1, T(1),
                        x+NBL1, 1);
        trsv_complex_up_ct(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i);

        IndexType j=0;

        for (; j+3<i; j+=4) {

            real_x0_.fill(real(x[j]));
            imag_x0_.fill(imag(x[j]));
            real_x1_.fill(real(x[j+1]));
            imag_x1_.fill(imag(x[j+1]));
            real_x2_.fill(real(x[j+2]));
            imag_x2_.fill(imag(x[j+2]));
            real_x3_.fill(real(x[j+3]));
            imag_x3_.fill(imag(x[j+3]));

            A0_.loadu(A+i+j*ldA);
            A1_.loadu(A+i+(j+1)*ldA);
            A2_.loadu(A+i+(j+2)*ldA);
            A3_.loadu(A+i+(j+3)*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j<i; ++j) {

            real_x0_.fill(real(x[j]));
            imag_x0_.fill(imag(x[j]));

            A0_.loadu(A+i+j*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, imag_x0_));
        }

        x_.storeu(x+i);


        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=i; j<i+ii; ++j) {
                x[i+ii] -= conj(A[(i+ii)+j*ldA])*x[j];
            }

            if (diag==NonUnit) {
                x[i+ii] /= conj(A[(i+ii)*(ldA+1)]);
            }
        }
    }

    for (; i<n; ++i) {

        for (IndexType j=0; j<i; ++j) {
            x[i] -= conj(A[i+j*ldA])*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= conj(A[i*(ldA+1)]);
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_lo_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_lo_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_complex_lo_n(diag, NBL3, A, ldA, x, 1);
        gemv_complex_n(n-NBL3, NBL3, T(-1), A+NBL3*ldA, ldA, x, 1, T(1),
                       x+NBL3, 1);
        trsv_complex_lo_n(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;
    } else if (n>NBL2 ) {

        trsv_complex_lo_n(diag, NBL2, A, ldA, x, 1);
        gemv_complex_n(n-NBL2, NBL2, T(-1), A+NBL2*ldA, ldA, x, 1, T(1),
                       x+NBL2, 1);
        trsv_complex_lo_n(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;
    } else if (n>NBL1 ) {

        trsv_complex_lo_n(diag, NBL1, A, ldA, x, 1);
        gemv_complex_n(n-NBL1, NBL1, T(-1), A+NBL1*ldA, ldA, x, 1, T(1),
                       x+NBL1, 1);
        trsv_complex_lo_n(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=0; i<n; ++i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        IntrinsicType result_;
        result_.setZero();

        T result[numElements];

        for (; j+4*numElements-1<i; j+=4*numElements) {

            A0_.loadu(A+i*ldA+j);
            A1_.loadu(A+i*ldA+j+numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

            x0_.loadu(x+j);
            x1_.loadu(x+j+numElements);
            x2_.loadu(x+j+2*numElements);
            x3_.loadu(x+j+3*numElements);

            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);
            real_x1_ = intrinsic_real_(x1_);
            imag_x1_ = intrinsic_imag_(x1_);
            real_x2_ = intrinsic_real_(x2_);
            imag_x2_ = intrinsic_imag_(x2_);
            real_x3_ = intrinsic_real_(x3_);
            imag_x3_ = intrinsic_imag_(x3_);


            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, real_x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, real_x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A1_, imag_x1_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A2_, imag_x2_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j+numElements-1<i; j+=numElements) {
            A0_.loadu(A+i*ldA+j);

            x0_.loadu(x+j);
            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }
        for (; j<i; ++j) {
            x[i] -= A[i*ldA+j]*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_lo_c(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_real_complex_lo_c");

    using std::conj;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_complex_lo_c(diag, NBL3, A, ldA, x, 1);
        gemv_complex_c(n-NBL3, NBL3, T(-1), A+NBL3*ldA, ldA, x, 1, T(1),
                       x+NBL3, 1);
        trsv_complex_lo_c(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);

        return;
    } else if (n>NBL2 ) {

        trsv_complex_lo_c(diag, NBL2, A, ldA, x, 1);
        gemv_complex_c(n-NBL2, NBL2, T(-1), A+NBL2*ldA, ldA, x, 1, T(1),
                       x+NBL2, 1);
        trsv_complex_lo_c(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);

        return;
    } else if (n>NBL1 ) {

        trsv_complex_lo_c(diag, NBL1, A, ldA, x, 1);
        gemv_complex_c(n-NBL1, NBL1, T(-1), A+NBL1*ldA, ldA, x, 1, T(1),
                       x+NBL1, 1);
        trsv_complex_lo_c(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);

        return;
    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    for (IndexType i=0; i<n; ++i) {

        IndexType j=0;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        IntrinsicType result_;
        result_.setZero();

        IntrinsicPrimitiveType minusOne_;
        PT minusOne(-1);
        minusOne_.fill(minusOne);

        T result[numElements];

        for (; j+4*numElements-1<i; j+=4*numElements) {

            A0_.loadu(A+i*ldA+j);
            A1_.loadu(A+i*ldA+j+numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

            x0_.loadu(x+j);
            x1_.loadu(x+j+numElements);
            x2_.loadu(x+j+2*numElements);
            x3_.loadu(x+j+3*numElements);

            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);
            real_x1_ = intrinsic_mul_(minusOne_, intrinsic_real_(x1_));
            imag_x1_ = intrinsic_imag_(x1_);
            real_x2_ = intrinsic_mul_(minusOne_, intrinsic_real_(x2_));
            imag_x2_ = intrinsic_imag_(x2_);
            real_x3_ = intrinsic_mul_(minusOne_, intrinsic_real_(x3_));
            imag_x3_ = intrinsic_imag_(x3_);


            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A1_, real_x1_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A2_, real_x2_));
            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, imag_x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, imag_x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j+numElements-1<i; j+=numElements) {
            A0_.loadu(A+i*ldA+j);

            x0_.loadu(x+j);
            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] -= result[k];
        }
        for (; j<i; ++j) {
            x[i] -= conj(A[i*ldA+j])*x[j];
        }

        if (diag==NonUnit) {
            x[i] /= conj(A[i*(ldA+1)]);
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_lo_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_lo_t");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {
        trsv_complex_lo_t(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_complex_t(NBL3, n-NBL3, T(-1), A+(n-NBL3)*ldA, ldA,
                       x+n-NBL3, 1, T(1), x, 1);
        trsv_complex_lo_t(diag, n-NBL3, A, ldA, x, 1);
        return;
    } else if (n>NBL2 ) {
        trsv_complex_lo_t(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_complex_t(NBL2, n-NBL2, T(-1), A+(n-NBL2)*ldA, ldA,
                       x+n-NBL2, 1, T(1), x, 1);
        trsv_complex_lo_t(diag, n-NBL2, A, ldA, x, 1);
        return;
    } else if (n>NBL1 ) {
        trsv_complex_lo_t(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_complex_t(NBL1, n-NBL1, T(-1), A+(n-NBL1)*ldA, ldA,
                       x+n-NBL1, 1, T(1), x, 1);
        trsv_complex_lo_t(diag, n-NBL1, A, ldA, x, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=n-1;

    for (; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i-numElements+1);

        IndexType j=0;

        for (; j+3<n-i-1; j+=4) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);
            A1_.loadu(A+(i+j+2)*ldA+i-numElements+1);
            A2_.loadu(A+(i+j+3)*ldA+i-numElements+1);
            A3_.loadu(A+(i+j+4)*ldA+i-numElements+1);

            real_x0_.fill( real(x[i+j+1]));
            imag_x0_.fill(-imag(x[i+j+1]));
            real_x1_.fill( real(x[i+j+2]));
            imag_x1_.fill(-imag(x[i+j+2]));
            real_x2_.fill( real(x[i+j+3]));
            imag_x2_.fill(-imag(x[i+j+3]));
            real_x3_.fill( real(x[i+j+4]));
            imag_x3_.fill(-imag(x[i+j+4]));

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j<n-i-1; ++j) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);

            real_x0_.fill( real(x[i+j+1]));
            imag_x0_.fill(-imag(x[i+j+1]));

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
        }
        x_.storeu(x+i-numElements+1);

        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=0; j<ii; ++j) {
                x[i-ii] -= A[(i-ii)*(ldA+1)+(j+1)*ldA]*x[(i-ii)+j+1];
            }

            if (diag==NonUnit) {
                x[i-ii] /= A[(i-ii)*(ldA+1)];
            }
        }
    }

    for (; i>=0; --i) {

        for (IndexType j=0; j<n-i-1; ++j) {
            x[i] -= A[i*(ldA+1)+(j+1)*ldA]*x[i+j+1];
        }

        if (diag==NonUnit) {
            x[i] /= A[i*(ldA+1)];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trsv_complex_lo_ct(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_complex_lo_ct");


    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRSV, T, IndexType>::NBlockL1();

    if (n>NBL3 ) {

        trsv_complex_lo_ct(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA, x+n-NBL3, 1);
        gemv_complex_ct(NBL3, n-NBL3, T(-1), A+(n-NBL3)*ldA, ldA,
                        x+n-NBL3, 1, T(1), x, 1);
        trsv_complex_lo_ct(diag, n-NBL3, A, ldA, x, 1);
        return;

    } else if (n>NBL3 ) {

        trsv_complex_lo_ct(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA, x+n-NBL2, 1);
        gemv_complex_ct(NBL2, n-NBL2, T(-1), A+(n-NBL2)*ldA, ldA,
                        x+n-NBL2, 1, T(1), x, 1);
        trsv_complex_lo_ct(diag, n-NBL2, A, ldA, x, 1);
        return;

    } else if (n>NBL1 ) {

        trsv_complex_lo_ct(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA, x+n-NBL1, 1);
        gemv_complex_ct(NBL1, n-NBL1, T(-1), A+(n-NBL1)*ldA, ldA,
                        x+n-NBL1, 1, T(1), x, 1);
        trsv_complex_lo_ct(diag, n-NBL1, A, ldA, x, 1);
        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=n-1;

    for (; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i-numElements+1);

        IndexType j=0;

        for (; j+3<n-i-1; j+=4) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);
            A1_.loadu(A+(i+j+2)*ldA+i-numElements+1);
            A2_.loadu(A+(i+j+3)*ldA+i-numElements+1);
            A3_.loadu(A+(i+j+4)*ldA+i-numElements+1);

            real_x0_.fill( real(x[i+j+1]));
            imag_x0_.fill( imag(x[i+j+1]));
            real_x1_.fill( real(x[i+j+2]));
            imag_x1_.fill( imag(x[i+j+2]));
            real_x2_.fill( real(x[i+j+3]));
            imag_x2_.fill( imag(x[i+j+3]));
            real_x3_.fill( real(x[i+j+4]));
            imag_x3_.fill( imag(x[i+j+4]));

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_sub_(x_, intrinsic_mul_(A3_, imag_x3_));
        }

        for (; j<n-i-1; ++j) {

            A0_.loadu(A+(i+j+1)*ldA+i-numElements+1);

            real_x0_.fill( real(x[i+j+1]));
            imag_x0_.fill( imag(x[i+j+1]));

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_sub_(x_, intrinsic_mul_(A0_, imag_x0_));
        }
        x_.storeu(x+i-numElements+1);

        for (IndexType ii=0; ii<numElements; ++ii) {

            for (IndexType j=0; j<ii; ++j) {
                x[i-ii] -= conj(A[(i-ii)*(ldA+1)+(j+1)*ldA])*x[(i-ii)+j+1];
            }

            if (diag==NonUnit) {
                x[i-ii] /= conj(A[(i-ii)*(ldA+1)]);
            }
        }
    }

    for (; i>=0; --i) {

        for (IndexType j=0; j<n-i-1; ++j) {
            x[i] -= conj(A[i*(ldA+1)+(j+1)*ldA])*x[i+j+1];
        }

        if (diag==NonUnit) {
            x[i] /= conj(A[i*(ldA+1)]);
        }
    }
}


template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
trsv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const T *A, IndexType ldA,
     T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_intrinsics [real, " INTRINSIC_NAME "]");


    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        upLo = (upLo==Upper) ? Lower : Upper;
        trsv(RowMajor, upLo, transA, diag, n, A, ldA, x, incX);
        return;
    }

    if ( incX != 1) {

        cxxblas::trsv<IndexType, T, T>(RowMajor, Upper, transA, diag, n, A, ldA,
                                       x, incX);
        return;

    }

    if ( upLo==Upper) {

        if  (transA==NoTrans || transA==Conj) {

            trsv_real_up_n(diag, n, A, ldA, x, incX);

        } else if (transA==Trans || transA==ConjTrans) {

            trsv_real_up_t(diag, n, A, ldA, x, incX);

        }

    } else if (upLo==Lower) {

        if  (transA==NoTrans || transA==Conj) {

            trsv_real_lo_n(diag, n, A, ldA, x, incX);

        } else if (transA==Trans || transA==ConjTrans) {

            trsv_real_lo_t(diag, n, A, ldA, x, incX);

        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
trsv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const T *A, IndexType ldA,
     T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trsv_intrinsics [complex, " INTRINSIC_NAME "]");


    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        upLo = (upLo==Upper) ? Lower : Upper;
        trsv(RowMajor, upLo, transA, diag, n, A, ldA, x, incX);
        return;
    }

    if ( incX != 1) {

        cxxblas::trsv<IndexType, T, T>(RowMajor, Upper, transA, diag, n, A, ldA,
                                       x, incX);
        return;

    }

    if ( upLo==Upper) {

        if  ( transA==NoTrans ) {

            trsv_complex_up_n(diag, n, A, ldA, x, incX);

        } else if  ( transA==Conj ) {

            trsv_complex_up_c(diag, n, A, ldA, x, incX);

        } else if ( transA==Trans ) {

            trsv_complex_up_t(diag, n, A, ldA, x, incX);

        } else if ( transA==ConjTrans ) {

            trsv_complex_up_ct(diag, n, A, ldA, x, incX);

        }

    } else if (upLo==Lower) {

        if  ( transA==NoTrans ) {

            trsv_complex_lo_n(diag, n, A, ldA, x, incX);

        } else if  ( transA==Conj ) {

            trsv_complex_lo_c(diag, n, A, ldA, x, incX);

        } else if ( transA==Trans ) {

            trsv_complex_lo_t(diag, n, A, ldA, x, incX);

        } else if ( transA==ConjTrans ) {

            trsv_complex_lo_ct(diag, n, A, ldA, x, incX);

        }
    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRSV_TCC
