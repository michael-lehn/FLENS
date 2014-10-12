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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRMV_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRMV_TCC 1

#include <cxxstd/iostream.h>

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trmv_real_up_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_up_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n> NBL3) {
        trmv_real_up_n(diag, NBL3, A, ldA, x, 1);
        gemv_real_n(NBL3, n-NBL3, T(1), A+NBL3, ldA, x+NBL3, 1, T(1), x, 1);
        trmv_real_up_n(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;
    } else if ( n> NBL2) {
        trmv_real_up_n(diag, NBL2, A, ldA, x, 1);
        gemv_real_n(NBL2, n-NBL2, T(1), A+NBL2, ldA, x+NBL2, 1, T(1), x, 1);
        trmv_real_up_n(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;
    } else if ( n> NBL1) {
        trmv_real_up_n(diag, NBL1, A, ldA, x, 1);
        gemv_real_n(NBL1, n-NBL1, T(1), A+NBL1, ldA, x+NBL1, 1, T(1), x, 1);
        trmv_real_up_n(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;
    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=0; i<n; ++i) {

        if (diag==NonUnit)
            x[i] = A[i*ldA+i]*x[i];

        IndexType j=i+1;

        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType result_;

        result_.setZero();

        for (; j+4*numElements-1<n; j+=4*numElements) {

            x0_.loadu(x+j              );
            x1_.loadu(x+j+  numElements);
            x2_.loadu(x+j+2*numElements);
            x3_.loadu(x+j+3*numElements);

            A0_.loadu(A+i*ldA+j              );
            A1_.loadu(A+i*ldA+j+  numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, x1_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, x2_));
            result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, x3_));
        }

        for (; j+numElements-1<n; j+=numElements) {
           x0_.loadu(x+j);
           A0_.loadu(A+i*ldA+j);

           result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] += result[k];
        }

        for (; j<n; ++j) {
           x[i] += A[i*ldA+j]*x[j];
        }

    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trmv_real_up_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_up_t");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3 ) {

        trmv_real_up_t(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                       x+(n-NBL3), IndexType(1));
        gemv_real_t(n-NBL3, NBL3, T(1), A+(n-NBL3), ldA,
                    x, IndexType(1), T(1), x+(n-NBL3), IndexType(1));
        trmv_real_up_t(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL2 ) {

        trmv_real_up_t(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                       x+(n-NBL2), IndexType(1));
        gemv_real_t(n-NBL2, NBL2, T(1), A+(n-NBL2), ldA,
                    x, IndexType(1), T(1), x+(n-NBL2), IndexType(1));
        trmv_real_up_t(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1 ) {

        trmv_real_up_t(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                       x+(n-NBL1), IndexType(1));
        gemv_real_t(n-NBL1, NBL1, T(1), A+(n-NBL1), ldA,
                    x, IndexType(1), T(1), x+(n-NBL1), IndexType(1));
        trmv_real_up_t(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    IndexType i  =  n-1;

    for(; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;


        IndexType firstIndex = i;
        IndexType lastIndex  = i-numElements+1;

        for (IndexType ii=firstIndex; ii>=lastIndex; --ii) {

            if (diag==NonUnit) {
                x[ii] = A[ii+ii*ldA]*x[ii];
            }

            for (IndexType jj=lastIndex; jj<ii; ++jj) {
                x[ii] += A[ii+jj*ldA]*x[jj];
            }
        }

        x_.loadu(x+lastIndex*incX);

        IndexType ii=0;
        for (; ii+3<lastIndex; ii+=4) {

            x0_.fill(x[ii]);
            x1_.fill(x[ii+1]);
            x2_.fill(x[ii+2]);
            x3_.fill(x[ii+3]);

            A0_.loadu(A+lastIndex+ii*ldA);
            A1_.loadu(A+lastIndex+ii*ldA+ldA);
            A2_.loadu(A+lastIndex+ii*ldA+2*ldA);
            A3_.loadu(A+lastIndex+ii*ldA+3*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, x3_));
        }


        for (; ii<lastIndex; ++ii) {
            A0_.loadu(A+lastIndex+ii*ldA);
            x0_.fill(x[ii]);
            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, x0_));
        }

        x_.storeu(x+lastIndex);
    }

    for (; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = A[i+i*ldA]*x[i];
        }

        for (IndexType j=0; j<i; ++j) {
            x[i] += A[i+j*ldA]*x[j];
        }

    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trmv_real_lo_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_lo_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_real_lo_n(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                       x+(n-NBL3), IndexType(1));
        gemv_real_n(NBL3, n-NBL3, T(1), A+(n-NBL3)*ldA, ldA,
                    x, incX, T(1), x+(n-NBL3), IndexType(1));
        trmv_real_lo_n(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL2  ) {

        trmv_real_lo_n(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                       x+(n-NBL2), IndexType(1));
        gemv_real_n(NBL2, n-NBL2, T(1), A+(n-NBL2)*ldA, ldA,
                    x, incX, T(1), x+(n-NBL2), IndexType(1));
        trmv_real_lo_n(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1  ) {

        trmv_real_lo_n(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                       x+(n-NBL1), IndexType(1));
        gemv_real_n(NBL1, n-NBL1, T(1), A+(n-NBL1)*ldA, ldA,
                    x, incX, T(1), x+(n-NBL1), IndexType(1));
        trmv_real_lo_n(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=n-1; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = A[i*ldA+i]*x[i];
        }
        IndexType j=0;

        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType result_;
        result_.setZero();

        for (; j+4*numElements-1<i; j+=4*numElements) {

           x0_.loadu(x+j               );
           x1_.loadu(x+j+  numElements);
           x2_.loadu(x+j+2*numElements);
           x3_.loadu(x+j+3*numElements);

           A0_.loadu(A+i*ldA+j              );
           A1_.loadu(A+i*ldA+j+  numElements);
           A2_.loadu(A+i*ldA+j+2*numElements);
           A3_.loadu(A+i*ldA+j+3*numElements);

           result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
           result_ = intrinsic_add_(result_, intrinsic_mul_(A1_, x1_));
           result_ = intrinsic_add_(result_, intrinsic_mul_(A2_, x2_));
           result_ = intrinsic_add_(result_, intrinsic_mul_(A3_, x3_));
        }

        for (; j+numElements-1<i; j+=numElements) {
           x0_.loadu(x+j);
           A0_.loadu(A+i*ldA+j);
           result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, x0_));
        }

        result_.storeu(&result[0]);

        for(IndexType k=0; k<numElements; ++k) {
           x[i] += result[k];
        }

        for (; j<i; ++j) {
           x[i] += A[i*ldA+j]*x[j];
        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
trmv_real_lo_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_lo_t");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_real_lo_t(diag, NBL3, A, ldA, x, 1);
        gemv_real_t(n-NBL3, NBL3, T(1), A+NBL3*ldA, ldA, x+NBL3, 1, T(1), x, 1);
        trmv_real_lo_t(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;

    } else if ( n>NBL2  ) {

        trmv_real_lo_t(diag, NBL2, A, ldA, x, 1);
        gemv_real_t(n-NBL2, NBL2, T(1), A+NBL2*ldA, ldA, x+NBL2, 1, T(1), x, 1);
        trmv_real_lo_t(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;

    } else if ( n>NBL1  ) {

        trmv_real_lo_t(diag, NBL1, A, ldA, x, 1);
        gemv_real_t(n-NBL1, NBL1, T(1), A+NBL1*ldA, ldA, x+NBL1, 1, T(1), x, 1);
        trmv_real_lo_t(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;

    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IndexType ii=0;

        for (; ii<numElements; ++ii) {

            if (diag==NonUnit) {
                x[i+ii] = A[i+ii+(i+ii)*ldA]*x[i+ii];

            }

            for (IndexType j=ii+1;j<numElements; ++j) {
               x[i+ii] += A[i+ii+(j+i)*ldA]*x[i+j];
            }
        }

        IntrinsicType x_;

        IntrinsicType A0_, A1_, A2_, A3_;

        IntrinsicType x0_, x1_, x2_, x3_;

        x_.loadu(x+i);

        for (; ii+i+3<n; ii+=4) {

            x0_.fill(x[ii+i  ]);
            x1_.fill(x[ii+i+1]);
            x2_.fill(x[ii+i+2]);
            x3_.fill(x[ii+i+3]);

            A0_.loadu(A+i+(ii+i  )*ldA);
            A1_.loadu(A+i+(ii+i+1)*ldA);
            A2_.loadu(A+i+(ii+i+2)*ldA);
            A3_.loadu(A+i+(ii+i+3)*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, x3_));
        }


        for (; ii+i<n; ++ii) {

            A0_.loadu(A+i+(ii+i)*ldA);

            x0_.fill(x[ii+i]);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, x0_));

        }

        x_.storeu(x+i);

    }


    for (; i<n; ++i) {

        if (diag==NonUnit) {
            x[i] = A[i*ldA+i]*x[i];
        }

        for (IndexType j=1; j<n-i; ++j) {

            x[i] += A[i*ldA+i+j*ldA]*x[i+j];

        }

    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_up_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_up_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3 ) {

        trmv_complex_up_n(diag, NBL3, A, ldA, x, 1);
        gemv_complex_n(NBL3, n-NBL3, T(1), A+NBL3, ldA, x+NBL3, 1, T(1), x, 1);
        trmv_complex_up_n(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;

    } else if ( n>NBL2 ) {

        trmv_complex_up_n(diag, NBL2, A, ldA, x, 1);
        gemv_complex_n(NBL2, n-NBL2, T(1), A+NBL2, ldA, x+NBL2, 1, T(1), x, 1);
        trmv_complex_up_n(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;

    } else if ( n>NBL1 ) {

        trmv_complex_up_n(diag, NBL1, A, ldA, x, 1);
        gemv_complex_n(NBL1, n-NBL1, T(1), A+NBL1, ldA, x+NBL1, 1, T(1), x, 1);
        trmv_complex_up_n(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;

    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=0; i<n; ++i) {

        if (diag==NonUnit)
            x[i] = A[i*ldA+i]*x[i];

        IndexType j=i+1;

        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType result_;

        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        result_.setZero();

        for (; j+4*numElements-1<n; j+=4*numElements) {

            x0_.loadu(x+j              );
            x1_.loadu(x+j+  numElements);
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


            A0_.loadu(A+i*ldA+j              );
            A1_.loadu(A+i*ldA+j+  numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

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

        for (; j+numElements-1<n; j+=numElements) {

           x0_.loadu(x+j);

            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);

           A0_.loadu(A+i*ldA+j);

           result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] += result[k];
        }

        for (; j<n; ++j) {
           x[i] += A[i*ldA+j]*x[j];
        }

    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_up_c(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_complex_c");

    using std::conj;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_up_c(diag, NBL3, A, ldA, x, 1);
        gemv_complex_c(NBL3, n-NBL3, T(1), A+NBL3, ldA, x+NBL3, 1, T(1), x, 1);
        trmv_complex_up_c(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;

    } else if ( n>NBL2  ) {

        trmv_complex_up_c(diag, NBL2, A, ldA, x, 1);
        gemv_complex_c(NBL2, n-NBL2, T(1), A+NBL2, ldA, x+NBL2, 1, T(1), x, 1);
        trmv_complex_up_c(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;

    } else if ( n>NBL1  ) {

        trmv_complex_up_c(diag, NBL1, A, ldA, x, 1);
        gemv_complex_c(NBL1, n-NBL1, T(1), A+NBL1, ldA, x+NBL1, 1, T(1), x, 1);
        trmv_complex_up_c(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;

    }


    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=0; i<n; ++i) {

        if (diag==NonUnit)
            x[i] = conj(A[i*ldA+i])*x[i];

        IndexType j=i+1;

        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType result_;

        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;
        IntrinsicPrimitiveType minusOne_;

        PT minusOne(-1);
        minusOne_.fill(minusOne);

        result_.setZero();

        for (; j+4*numElements-1<n; j+=4*numElements) {

            x0_.loadu(x+j              );
            x1_.loadu(x+j+  numElements);
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


            A0_.loadu(A+i*ldA+j              );
            A1_.loadu(A+i*ldA+j+  numElements);
            A2_.loadu(A+i*ldA+j+2*numElements);
            A3_.loadu(A+i*ldA+j+3*numElements);

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

        for (; j+numElements-1<n; j+=numElements) {

           x0_.loadu(x+j);

            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);

           A0_.loadu(A+i*ldA+j);

           result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for (IndexType k=0; k<numElements; ++k) {
            x[i] += result[k];
        }

        for (; j<n; ++j) {
           x[i] += conj(A[i*ldA+j])*x[j];
        }

    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_up_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_up_t");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_up_t(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                          x+(n-NBL3), IndexType(1));
        gemv_complex_t(n-NBL3, NBL3, T(1), A+(n-NBL3), ldA,
                       x, IndexType(1), T(1), x+(n-NBL3), IndexType(1));
        trmv_complex_up_t(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL3  ) {

        trmv_complex_up_t(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                          x+(n-NBL2), IndexType(1));
        gemv_complex_t(n-NBL2, NBL2, T(1), A+(n-NBL2), ldA,
                       x, IndexType(1), T(1), x+(n-NBL2), IndexType(1));
        trmv_complex_up_t(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1  ) {

        trmv_complex_up_t(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                          x+(n-NBL1), IndexType(1));
        gemv_complex_t(n-NBL1, NBL1, T(1), A+(n-NBL1), ldA,
                       x, IndexType(1), T(1), x+(n-NBL1), IndexType(1));
        trmv_complex_up_t(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i  =  n-1;

    for(; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        IndexType firstIndex = i;
        IndexType lastIndex  = i-numElements+1;

        for (IndexType ii=firstIndex; ii>=lastIndex; --ii) {

            if (diag==NonUnit) {
                x[ii] = A[ii+ii*ldA]*x[ii];
            }

            for (IndexType jj=lastIndex; jj<ii; ++jj) {
                x[ii] += A[ii+jj*ldA]*x[jj];
            }
        }

        x_.loadu(x+lastIndex*incX);

        IndexType ii=0;
        for (; ii+3<lastIndex; ii+=4) {

            real_x0_.fill(real(x[ii]));
            real_x1_.fill(real(x[ii+1]));
            real_x2_.fill(real(x[ii+2]));
            real_x3_.fill(real(x[ii+3]));

            imag_x0_.fill(imag(x[ii]));
            imag_x1_.fill(imag(x[ii+1]));
            imag_x2_.fill(imag(x[ii+2]));
            imag_x3_.fill(imag(x[ii+3]));

            A0_.loadu(A+lastIndex+ii*ldA);
            A1_.loadu(A+lastIndex+ii*ldA+ldA);
            A2_.loadu(A+lastIndex+ii*ldA+2*ldA);
            A3_.loadu(A+lastIndex+ii*ldA+3*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, imag_x3_));
        }


        for (; ii<lastIndex; ++ii) {

            real_x0_.fill(real(x[ii]));
            imag_x0_.fill(imag(x[ii]));

            A0_.loadu(A+lastIndex+ii*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
        }

        x_.storeu(x+lastIndex);
    }

    for (; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = A[i+i*ldA]*x[i];
        }

        for (IndexType j=0; j<i; ++j) {
            x[i] += A[i+j*ldA]*x[j];
        }

    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_up_ct(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_up_ct");

    using std::conj;
    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_up_ct(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                           x+(n-NBL3), IndexType(1));
        gemv_complex_ct(n-NBL3, NBL3, T(1), A+(n-NBL3), ldA,
                        x, IndexType(1), T(1), x+(n-NBL3), IndexType(1));
        trmv_complex_up_ct(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL2  ) {

        trmv_complex_up_ct(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                           x+(n-NBL2), IndexType(1));
        gemv_complex_ct(n-NBL2, NBL2, T(1), A+(n-NBL2), ldA,
                        x, IndexType(1), T(1), x+(n-NBL2), IndexType(1));
        trmv_complex_up_ct(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1  ) {

        trmv_complex_up_ct(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                           x+(n-NBL1), IndexType(1));
        gemv_complex_ct(n-NBL1, NBL1, T(1), A+(n-NBL1), ldA,
                        x, IndexType(1), T(1), x+(n-NBL1), IndexType(1));
        trmv_complex_up_ct(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    IndexType i  =  n-1;

    for(; i-numElements+1>=0; i-=numElements) {

        IntrinsicType x_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        IndexType firstIndex = i;
        IndexType lastIndex  = i-numElements+1;

        for (IndexType ii=firstIndex; ii>=lastIndex; --ii) {

            if (diag==NonUnit) {
                x[ii] = conj(A[ii+ii*ldA])*x[ii];
            }

            for (IndexType jj=lastIndex; jj<ii; ++jj) {
                x[ii] += conj(A[ii+jj*ldA])*x[jj];
            }
        }

        x_.loadu(x+lastIndex*incX);

        IndexType ii=0;
        for (; ii+3<lastIndex; ii+=4) {

            real_x0_.fill(-real(x[ii]));
            real_x1_.fill(-real(x[ii+1]));
            real_x2_.fill(-real(x[ii+2]));
            real_x3_.fill(-real(x[ii+3]));

            imag_x0_.fill(imag(x[ii]));
            imag_x1_.fill(imag(x[ii+1]));
            imag_x2_.fill(imag(x[ii+2]));
            imag_x3_.fill(imag(x[ii+3]));

            A0_.loadu(A+lastIndex+ii*ldA);
            A1_.loadu(A+lastIndex+ii*ldA+ldA);
            A2_.loadu(A+lastIndex+ii*ldA+2*ldA);
            A3_.loadu(A+lastIndex+ii*ldA+3*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, imag_x3_));
        }


        for (; ii<lastIndex; ++ii) {

            real_x0_.fill(-real(x[ii]));
            imag_x0_.fill( imag(x[ii]));

            A0_.loadu(A+lastIndex+ii*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, imag_x0_));
        }

        x_.storeu(x+lastIndex);
    }

    for (; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = conj(A[i+i*ldA])*x[i];
        }

        for (IndexType j=0; j<i; ++j) {
            x[i] += conj(A[i+j*ldA])*x[j];
        }

    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_lo_n(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_lo_n");

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_lo_n(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                          x+(n-NBL3), IndexType(1));
        gemv_complex_n(NBL3, n-NBL3, T(1), A+(n-NBL3)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL3), IndexType(1));
        trmv_complex_lo_n(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL2  ) {

        trmv_complex_lo_n(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                          x+(n-NBL2), IndexType(1));
        gemv_complex_n(NBL2, n-NBL2, T(1), A+(n-NBL2)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL2), IndexType(1));
        trmv_complex_lo_n(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1  ) {

        trmv_complex_lo_n(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                          x+(n-NBL1), IndexType(1));
        gemv_complex_n(NBL1, n-NBL1, T(1), A+(n-NBL1)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL1), IndexType(1));
        trmv_complex_lo_n(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=n-1; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = A[i*ldA+i]*x[i];
        }
        IndexType j=0;

        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;
        IntrinsicType result_;
        result_.setZero();

        for (; j+4*numElements-1<i; j+=4*numElements) {

           x0_.loadu(x+j               );
           x1_.loadu(x+j+  numElements);
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

           A0_.loadu(A+i*ldA+j              );
           A1_.loadu(A+i*ldA+j+  numElements);
           A2_.loadu(A+i*ldA+j+2*numElements);
           A3_.loadu(A+i*ldA+j+3*numElements);

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

           x0_.loadu(x+j);

            real_x0_ = intrinsic_real_(x0_);
            imag_x0_ = intrinsic_imag_(x0_);

           A0_.loadu(A+i*ldA+j);

           result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for(IndexType k=0; k<numElements; ++k) {
           x[i] += result[k];
        }

        for (; j<i; ++j) {
           x[i] += A[i*ldA+j]*x[j];
        }
    }

}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_lo_c(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_real_complex_lo_c");

    using std::conj;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_lo_c(diag, NBL3, A+(n-NBL3)*(ldA+1), ldA,
                          x+(n-NBL3), IndexType(1));
        gemv_complex_c(NBL3, n-NBL3, T(1), A+(n-NBL3)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL3), IndexType(1));
        trmv_complex_lo_c(diag, n-NBL3, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL2  ) {

        trmv_complex_lo_c(diag, NBL2, A+(n-NBL2)*(ldA+1), ldA,
                          x+(n-NBL2), IndexType(1));
        gemv_complex_c(NBL2, n-NBL2, T(1), A+(n-NBL2)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL2), IndexType(1));
        trmv_complex_lo_c(diag, n-NBL2, A, ldA, x, IndexType(1));

        return;

    } else if ( n>NBL1  ) {

        trmv_complex_lo_c(diag, NBL1, A+(n-NBL1)*(ldA+1), ldA,
                          x+(n-NBL1), IndexType(1));
        gemv_complex_c(NBL1, n-NBL1, T(1), A+(n-NBL1)*ldA, ldA,
                       x, incX, T(1), x+(n-NBL1), IndexType(1));
        trmv_complex_lo_c(diag, n-NBL1, A, ldA, x, IndexType(1));

        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    T result[numElements];

    for (IndexType i=n-1; i>=0; --i) {

        if (diag==NonUnit) {
            x[i] = conj(A[i*ldA+i])*x[i];
        }
        IndexType j=0;

        IntrinsicType x0_, x1_, x2_, x3_;
        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;
        IntrinsicType result_;


        result_.setZero();

        IntrinsicPrimitiveType minusOne_;
        PT minusOne(-1);
        minusOne_.fill(minusOne);


        for (; j+4*numElements-1<i; j+=4*numElements) {

           x0_.loadu(x+j               );
           x1_.loadu(x+j+  numElements);
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

           A0_.loadu(A+i*ldA+j              );
           A1_.loadu(A+i*ldA+j+  numElements);
           A2_.loadu(A+i*ldA+j+2*numElements);
           A3_.loadu(A+i*ldA+j+3*numElements);

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

            x0_.loadu(x+j);

            real_x0_ = intrinsic_mul_(minusOne_, intrinsic_real_(x0_));
            imag_x0_ = intrinsic_imag_(x0_);

            A0_.loadu(A+i*ldA+j);

            result_ = intrinsic_addsub_(result_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            result_ = intrinsic_add_(result_, intrinsic_mul_(A0_, imag_x0_));
        }

        result_.storeu(&result[0]);

        for(IndexType k=0; k<numElements; ++k) {
           x[i] += result[k];
        }

        for (; j<i; ++j) {
           x[i] += conj(A[i*ldA+j])*x[j];
        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_lo_t(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_lo_t");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_lo_t(diag, NBL3, A, ldA, x, 1);
        gemv_complex_t(n-NBL3, NBL3, T(1), A+NBL3*ldA, ldA, x+NBL3, 1, T(1),
                       x, 1);
        trmv_complex_lo_t(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;

    } else if ( n>NBL2  ) {

        trmv_complex_lo_t(diag, NBL2, A, ldA, x, 1);
        gemv_complex_t(n-NBL2, NBL2, T(1), A+NBL2*ldA, ldA, x+NBL2, 1, T(1),
                       x, 1);
        trmv_complex_lo_t(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;

    } else if ( n>NBL1  ) {

        trmv_complex_lo_t(diag, NBL1, A, ldA, x, 1);
        gemv_complex_t(n-NBL1, NBL1, T(1), A+NBL1*ldA, ldA, x+NBL1, 1, T(1),
                       x, 1);
        trmv_complex_lo_t(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;


    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IndexType ii=0;

        for (; ii<numElements; ++ii) {

            if (diag==NonUnit) {
                x[i+ii] = A[i+ii+(i+ii)*ldA]*x[i+ii];

            }

            for (IndexType j=ii+1;j<numElements; ++j) {
               x[i+ii] += A[i+ii+(j+i)*ldA]*x[i+j];
            }
        }

        IntrinsicType x_;

        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i);

        for (; ii+i+3<n; ii+=4) {

            real_x0_.fill(real(x[ii+i  ]));
            real_x1_.fill(real(x[ii+i+1]));
            real_x2_.fill(real(x[ii+i+2]));
            real_x3_.fill(real(x[ii+i+3]));

            imag_x0_.fill(imag(x[ii+i  ]));
            imag_x1_.fill(imag(x[ii+i+1]));
            imag_x2_.fill(imag(x[ii+i+2]));
            imag_x3_.fill(imag(x[ii+i+3]));

            A0_.loadu(A+i+(ii+i  )*ldA);
            A1_.loadu(A+i+(ii+i+1)*ldA);
            A2_.loadu(A+i+(ii+i+2)*ldA);
            A3_.loadu(A+i+(ii+i+3)*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, imag_x3_));

        }


        for (; ii+i<n; ++ii) {

            real_x0_.fill(real(x[ii+i]));
            imag_x0_.fill(imag(x[ii+i]));

            A0_.loadu(A+i+(ii+i)*ldA);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, imag_x0_));

        }

        x_.storeu(x+i);

    }


    for (; i<n; ++i) {

        if (diag==NonUnit) {
            x[i] = A[i*ldA+i]*x[i];
        }

        for (IndexType j=1; j<n-i; ++j) {

            x[i] += A[i*ldA+i+j*ldA]*x[i+j];

        }

    }


}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
trmv_complex_lo_ct(Diag diag, IndexType n,
               const T *A, IndexType ldA, T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_complex_lo_ct");


    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    const IndexType NBL3 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL3();
    const IndexType NBL2 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL2();
    const IndexType NBL1 = BlockSize<BlasFct::TRMV, T, IndexType>::NBlockL1();

    if ( n>NBL3  ) {

        trmv_complex_lo_ct(diag, NBL3, A, ldA, x, 1);
        gemv_complex_ct(n-NBL3, NBL3, T(1), A+NBL3*ldA, ldA, x+NBL3, 1, T(1),
                        x, 1);
        trmv_complex_lo_ct(diag, n-NBL3, A+NBL3*(ldA+1), ldA, x+NBL3, 1);
        return;

    } else if ( n>NBL2  ) {

        trmv_complex_lo_ct(diag, NBL2, A, ldA, x, 1);
        gemv_complex_ct(n-NBL2, NBL2, T(1), A+NBL2*ldA, ldA, x+NBL2, 1, T(1),
                        x, 1);
        trmv_complex_lo_ct(diag, n-NBL2, A+NBL2*(ldA+1), ldA, x+NBL2, 1);
        return;

    } else if ( n>NBL1  ) {

        trmv_complex_lo_ct(diag, NBL1, A, ldA, x, 1);
        gemv_complex_ct(n-NBL1, NBL1, T(1), A+NBL1*ldA, ldA, x+NBL1, 1, T(1),
                        x, 1);
        trmv_complex_lo_ct(diag, n-NBL1, A+NBL1*(ldA+1), ldA, x+NBL1, 1);
        return;

    }

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;


    IndexType i=0;

    for (; i+numElements-1<n; i+=numElements) {

        IndexType ii=0;

        for (; ii<numElements; ++ii) {

            if (diag==NonUnit) {
                x[i+ii] = conj(A[i+ii+(i+ii)*ldA])*x[i+ii];

            }

            for (IndexType j=ii+1;j<numElements; ++j) {
               x[i+ii] += conj(A[i+ii+(j+i)*ldA])*x[i+j];
            }
        }

        IntrinsicType x_;

        IntrinsicType A0_, A1_, A2_, A3_;
        IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
        IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

        x_.loadu(x+i);

        for (; ii+i+3<n; ii+=4) {

            real_x0_.fill(-real(x[ii+i  ]));
            real_x1_.fill(-real(x[ii+i+1]));
            real_x2_.fill(-real(x[ii+i+2]));
            real_x3_.fill(-real(x[ii+i+3]));

            imag_x0_.fill( imag(x[ii+i  ]));
            imag_x1_.fill( imag(x[ii+i+1]));
            imag_x2_.fill( imag(x[ii+i+2]));
            imag_x3_.fill( imag(x[ii+i+3]));

            A0_.loadu(A+i+(ii+i  )*ldA);
            A1_.loadu(A+i+(ii+i+1)*ldA);
            A2_.loadu(A+i+(ii+i+2)*ldA);
            A3_.loadu(A+i+(ii+i+3)*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A1_, real_x1_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A2_, real_x2_));
            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A3_, real_x3_));

            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, imag_x0_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A1_, imag_x1_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A2_, imag_x2_));
            x_ = intrinsic_add_(x_, intrinsic_mul_(A3_, imag_x3_));

        }


        for (; ii+i<n; ++ii) {

            real_x0_.fill(-real(x[ii+i]));
            imag_x0_.fill( imag(x[ii+i]));

            A0_.loadu(A+i+(ii+i)*ldA);

            x_ = intrinsic_addsub_(x_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            x_ = intrinsic_add_(x_, intrinsic_mul_(A0_, imag_x0_));

        }

        x_.storeu(x+i);

    }


    for (; i<n; ++i) {

        if (diag==NonUnit) {
            x[i] = conj(A[i*ldA+i])*x[i];
        }

        for (IndexType j=1; j<n-i; ++j) {

            x[i] += conj(A[i*ldA+i+j*ldA])*x[i+j];

        }

    }

}


template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const T *A, IndexType ldA,
     T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_intrinsics [real, " INTRINSIC_NAME "]");


    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        upLo = (upLo==Upper) ? Lower : Upper;
        trmv(RowMajor, upLo, transA, diag, n, A, ldA, x, incX);
        return;
    }

    if ( incX != 1) {

        cxxblas::trmv<IndexType, T, T>(RowMajor, Upper, transA, diag, n,
                                       A, ldA, x, incX);
        return;

    }

    if ( upLo==Upper) {

        if  (transA==NoTrans || transA==Conj) {

            trmv_real_up_n(diag, n, A, ldA, x, incX);

        } else if (transA==Trans || transA==ConjTrans) {

            trmv_real_up_t(diag, n, A, ldA, x, incX);

        }

    } else if (upLo==Lower) {

        if  (transA==NoTrans || transA==Conj) {

            trmv_real_lo_n(diag, n, A, ldA, x, incX);

        } else if (transA==Trans || transA==ConjTrans) {

            trmv_real_lo_t(diag, n, A, ldA, x, incX);

        }
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
trmv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const T *A, IndexType ldA,
     T *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("trmv_intrinsics [complex, " INTRINSIC_NAME "]");


    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        upLo = (upLo==Upper) ? Lower : Upper;
        trmv(RowMajor, upLo, transA, diag, n, A, ldA, x, incX);
        return;
    }

    if ( incX != 1) {

        cxxblas::trmv<IndexType, T, T>(RowMajor, Upper, transA, diag, n,
                                       A, ldA, x, incX);
        return;

    }

    if ( upLo==Upper) {

        if  ( transA==NoTrans ) {

            trmv_complex_up_n(diag, n, A, ldA, x, incX);

        } else if  ( transA==Conj ) {

            trmv_complex_up_c(diag, n, A, ldA, x, incX);

        } else if ( transA==Trans ) {

            trmv_complex_up_t(diag, n, A, ldA, x, incX);

        } else if ( transA==ConjTrans ) {

            trmv_complex_up_ct(diag, n, A, ldA, x, incX);

        }

    } else if (upLo==Lower) {

        if  ( transA==NoTrans ) {

            trmv_complex_lo_n(diag, n, A, ldA, x, incX);

        } else if  ( transA==Conj ) {

            trmv_complex_lo_c(diag, n, A, ldA, x, incX);

        } else if ( transA==Trans ) {

            trmv_complex_lo_t(diag, n, A, ldA, x, incX);

        } else if ( transA==ConjTrans ) {

            trmv_complex_lo_ct(diag, n, A, ldA, x, incX);

        }
    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_TRMV_TCC
