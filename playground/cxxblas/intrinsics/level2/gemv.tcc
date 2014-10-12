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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_GEMV_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_GEMV_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    gemv_real_n(IndexType m, IndexType n,
                const T &alpha,
                const T *A, IndexType ldA,
                const T *x, IndexType incX,
                const T &beta,
                T *y, IndexType incY)

{
    CXXBLAS_DEBUG_OUT("gemv_real_n");

    ASSERT( incX==1 );

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;

    scal(m, beta, y, incY);

    if (incY<0) {
        y -= incY*(m-1);
    }

    IndexType i=0;
    IndexType iY=0;

    T tmp_result[4][numElements];

     for (; i+3<m; i+=4, iY+=4*incY) {

        IntrinsicType A0_, A1_, A2_, A3_, A4_;
        IntrinsicType x_;
        IntrinsicType y0_, y1_, y2_, y3_;

        y0_.setZero();
        y1_.setZero();
        y2_.setZero();
        y3_.setZero();


        IndexType j=0;

        for (; j+numElements-1<n; j+=numElements) {
            A0_.loadu(A+i*ldA+j      );
            A1_.loadu(A+i*ldA+j+  ldA);
            A2_.loadu(A+i*ldA+j+2*ldA);
            A3_.loadu(A+i*ldA+j+3*ldA);

            x_.loadu(x+j);
            y0_ = intrinsic_add_(y0_, intrinsic_mul_(A0_, x_));
            y1_ = intrinsic_add_(y1_, intrinsic_mul_(A1_, x_));
            y2_ = intrinsic_add_(y2_, intrinsic_mul_(A2_, x_));
            y3_ = intrinsic_add_(y3_, intrinsic_mul_(A3_, x_));

        }

        y0_.storeu(&tmp_result[0][0]);
        y1_.storeu(&tmp_result[1][0]);
        y2_.storeu(&tmp_result[2][0]);
        y3_.storeu(&tmp_result[3][0]);

        for (;j<n;++j) {
            tmp_result[0][0] += A[i*ldA+j      ]*x[j];
            tmp_result[1][0] += A[i*ldA+j+  ldA]*x[j];
            tmp_result[2][0] += A[i*ldA+j+2*ldA]*x[j];
            tmp_result[3][0] += A[i*ldA+j+3*ldA]*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0][0] += tmp_result[0][k+1];
            tmp_result[1][0] += tmp_result[1][k+1];
            tmp_result[2][0] += tmp_result[2][k+1];
            tmp_result[3][0] += tmp_result[3][k+1];
        }

        y[iY       ] += alpha*tmp_result[0][0];
        y[iY+  incY] += alpha*tmp_result[1][0];
        y[iY+2*incY] += alpha*tmp_result[2][0];
        y[iY+3*incY] += alpha*tmp_result[3][0];
    }


    for (; i<m; ++i, iY+=incY) {

        IntrinsicType A_, x_, y_;
        y_.setZero();

        IndexType j=0;

        for (; j+numElements-1<n; j+=numElements) {
            A_.loadu(A+i*ldA+j);
            x_.loadu(x+j);
            y_ = intrinsic_add_(y_, intrinsic_mul_(A_, x_));
        }
        T tmp_result[numElements];
        y_.storeu(tmp_result);
        for (;j<n;++j) {
            tmp_result[0] += A[i*ldA+j]*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0] += tmp_result[k+1];
        }
        y[iY] += alpha*tmp_result[0];
    }

}

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
    gemv_real_t(IndexType m, IndexType n,
                const T &alpha,
                const T *A, IndexType ldA,
                const T *x, IndexType incX,
                const T &beta,
                T *y, IndexType incY)

{

    CXXBLAS_DEBUG_OUT("gemv_real_t");

    ASSERT( incY==1);

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
    const int numElements = IntrinsicType::numElements;


    IntrinsicType y_;
    IntrinsicType beta_(beta);

    IntrinsicType A0_, A1_, A2_, A3_;

    IntrinsicType x0_, x1_, x2_, x3_;

    scal(n, beta, y, 1);

    if (incX<0) {
        x -= incX*(m-1);
    }

    IndexType j=0, jX=0;

    for (; j+3<m; j+=4, jX+=4*incX) {

        x0_.fill(alpha*x[jX       ]);
        x1_.fill(alpha*x[jX+1*incX]);
        x2_.fill(alpha*x[jX+2*incX]);
        x3_.fill(alpha*x[jX+3*incX]);

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );
            A1_.loadu(A+i+j*ldA+  ldA);
            A2_.loadu(A+i+j*ldA+2*ldA);
            A3_.loadu(A+i+j*ldA+3*ldA);

            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, x0_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A1_, x1_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A2_, x2_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A3_, x3_));

            y_.storeu(y+i);
        }
    }

    for (; j<m; ++j, jX+=incX) {

        x0_.fill(alpha*x[jX]);

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );

            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, x0_));

            y_.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);

    for (; i<n; ++i) {
        T y_;
        dot(m, A+i, ldA, x, incX, y_);
        y[i] += alpha*y_;
    }

}

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
    gemv_complex_n(IndexType m, IndexType n,
                   const T &alpha,
                   const T *A, IndexType ldA,
                   const T *x, IndexType incX,
                   const T &beta,
                   T *y, IndexType incY)

{
    CXXBLAS_DEBUG_OUT("gemv_complex_n");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    scal(m, beta, y, incY);

    if (incY<0) {
        y -= incY*(m-1);
    }

    IntrinsicType A0_, A1_, A2_, A3_;
    IntrinsicType  x_;
    IntrinsicType  y0_, y1_, y2_, y3_;
    IntrinsicPrimitiveType real_x_, imag_x_;

    T tmp_result[8][numElements];

    IndexType i=0, iY=0;
    for (; i+3<m; i+=4, iY+=4*incY) {


    y0_.setZero();
    y1_.setZero();
    y2_.setZero();
    y3_.setZero();

    IndexType j=0;
    for (; j+numElements-1<n; j+=numElements) {

        A0_.loadu(A+i*ldA+j);
        A1_.loadu(A+i*ldA+j+ldA);
        A2_.loadu(A+i*ldA+j+2*ldA);
        A3_.loadu(A+i*ldA+j+3*ldA);

        x_.loadu(x+j);

        real_x_ = intrinsic_real_(x_);
        imag_x_ = intrinsic_imag_(x_);

        y0_ = intrinsic_add_(y0_, intrinsic_mul_(A0_, real_x_));
        y1_ = intrinsic_add_(y1_, intrinsic_mul_(A1_, real_x_));
        y2_ = intrinsic_add_(y2_, intrinsic_mul_(A2_, real_x_));
        y3_ = intrinsic_add_(y3_, intrinsic_mul_(A3_, real_x_));


        A0_ = intrinsic_swap_real_imag_(A0_);
        A1_ = intrinsic_swap_real_imag_(A1_);
        A2_ = intrinsic_swap_real_imag_(A2_);
        A3_ = intrinsic_swap_real_imag_(A3_);

        y0_ = intrinsic_addsub_(y0_, intrinsic_mul_(A0_, imag_x_));
        y1_ = intrinsic_addsub_(y1_, intrinsic_mul_(A1_, imag_x_));
        y2_ = intrinsic_addsub_(y2_, intrinsic_mul_(A2_, imag_x_));
        y3_ = intrinsic_addsub_(y3_, intrinsic_mul_(A3_, imag_x_));

        }

        y0_.storeu(&tmp_result[0][0]);
        y1_.storeu(&tmp_result[1][0]);
        y2_.storeu(&tmp_result[2][0]);
        y3_.storeu(&tmp_result[3][0]);

        for (;j<n;++j) {
            tmp_result[0][0] += A[i*ldA+j      ]*x[j];
            tmp_result[1][0] += A[i*ldA+j+  ldA]*x[j];
            tmp_result[2][0] += A[i*ldA+j+2*ldA]*x[j];
            tmp_result[3][0] += A[i*ldA+j+3*ldA]*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0][0] += tmp_result[0][k+1];
            tmp_result[1][0] += tmp_result[1][k+1];
            tmp_result[2][0] += tmp_result[2][k+1];
            tmp_result[3][0] += tmp_result[3][k+1];
        }

        y[iY       ] += alpha*tmp_result[0][0];
        y[iY+1*incY] += alpha*tmp_result[1][0];
        y[iY+2*incY] += alpha*tmp_result[2][0];
        y[iY+3*incY] += alpha*tmp_result[3][0];

    }



    for (; i<m; ++i, iY+=incY) {

        y0_.setZero();
        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            A0_.loadu(A+i*ldA+j);
            x_.loadu(x+j);

            real_x_ = intrinsic_real_(x_);
            imag_x_ = intrinsic_imag_(x_);


            y0_ = intrinsic_add_(y0_, intrinsic_mul_(A0_, real_x_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            y0_ = intrinsic_addsub_(y0_, intrinsic_mul_(A0_, imag_x_));
        }

        y0_.storeu(&tmp_result[0][0]);

        for (;j<n;++j) {
            tmp_result[0][0] += A[i*ldA+j    ]*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0][0] += tmp_result[0][k+1];
        }
        y[iY] += alpha*tmp_result[0][0];
    }


}

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
    gemv_complex_c(IndexType m, IndexType n,
                   const T &alpha,
                   const T *A, IndexType ldA,
                   const T *x, IndexType incX,
                   const T &beta,
                   T *y, IndexType incY)

{
    CXXBLAS_DEBUG_OUT("gemv_complex_c");

    using std::imag;
    using std::real;

    ASSERT( incX==1 );

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    scal(m, beta, y, incY);

    if (incY<0) {
        y -= incY*(m-1);
    }

    IntrinsicType A0_, A1_, A2_, A3_;
    IntrinsicType  x_;
    IntrinsicType  y0_, y1_, y2_, y3_;
    IntrinsicPrimitiveType real_x_, imag_x_, minusOne_;

    PT minusOne(-1);
    minusOne_.fill(minusOne);

    T tmp_result[4][numElements];

    IndexType i=0, iY=0;
    for (; i+3<m; i+=4, iY+=4*incY) {


        y0_.setZero();
        y1_.setZero();
        y2_.setZero();
        y3_.setZero();

        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            A0_.loadu(A+i*ldA+j);
            A1_.loadu(A+i*ldA+j+ldA);
            A2_.loadu(A+i*ldA+j+2*ldA);
            A3_.loadu(A+i*ldA+j+3*ldA);

            x_.loadu(x+j);

            real_x_ = intrinsic_mul_(minusOne_, intrinsic_real_(x_));
            imag_x_ = intrinsic_imag_(x_);

            y0_ = intrinsic_addsub_(y0_, intrinsic_mul_(A0_, real_x_));
            y1_ = intrinsic_addsub_(y1_, intrinsic_mul_(A1_, real_x_));
            y2_ = intrinsic_addsub_(y2_, intrinsic_mul_(A2_, real_x_));
            y3_ = intrinsic_addsub_(y3_, intrinsic_mul_(A3_, real_x_));


            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);

            y0_ = intrinsic_add_(y0_, intrinsic_mul_(A0_, imag_x_));
            y1_ = intrinsic_add_(y1_, intrinsic_mul_(A1_, imag_x_));
            y2_ = intrinsic_add_(y2_, intrinsic_mul_(A2_, imag_x_));
            y3_ = intrinsic_add_(y3_, intrinsic_mul_(A3_, imag_x_));

        }

        y0_.storeu(&tmp_result[0][0]);
        y1_.storeu(&tmp_result[1][0]);
        y2_.storeu(&tmp_result[2][0]);
        y3_.storeu(&tmp_result[3][0]);

        for (;j<n;++j) {
            tmp_result[0][0] += conj(A[i*ldA+j      ])*x[j];
            tmp_result[1][0] += conj(A[i*ldA+j+  ldA])*x[j];
            tmp_result[2][0] += conj(A[i*ldA+j+2*ldA])*x[j];
            tmp_result[3][0] += conj(A[i*ldA+j+3*ldA])*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0][0] += tmp_result[0][k+1];
            tmp_result[1][0] += tmp_result[1][k+1];
            tmp_result[2][0] += tmp_result[2][k+1];
            tmp_result[3][0] += tmp_result[3][k+1];
        }

        y[iY       ] += alpha*tmp_result[0][0];
        y[iY+1*incY] += alpha*tmp_result[1][0];
        y[iY+2*incY] += alpha*tmp_result[2][0];
        y[iY+3*incY] += alpha*tmp_result[3][0];

    }



    for (; i<m; ++i, iY+=incY) {

        y0_.setZero();
        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            A0_.loadu(A+i*ldA+j);
            x_.loadu(x+j);

            real_x_ = intrinsic_mul_(minusOne_, intrinsic_real_(x_));
            imag_x_ = intrinsic_imag_(x_);


            y0_ = intrinsic_addsub_(y0_, intrinsic_mul_(A0_, real_x_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            y0_ = intrinsic_add_(y0_, intrinsic_mul_(A0_, imag_x_));
        }

        y0_.storeu(&tmp_result[0][0]);

        for (;j<n;++j) {
            tmp_result[0][0] += conj(A[i*ldA+j    ])*x[j];
        }
        for (IndexType k=0; k+1<numElements; ++k) {
            tmp_result[0][0] += tmp_result[0][k+1];
        }
        y[iY] += alpha*tmp_result[0][0];
    }


}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
gemv_complex_t(IndexType m, IndexType n,
               const T &alpha,
               const T *A, IndexType ldA,
               const T *x, IndexType incX,
               const T &beta,
               T *y, IndexType incY)

{
    CXXBLAS_DEBUG_OUT("gemv_complex_t");

    using std::imag;
    using std::real;

    ASSERT( incY==1 );

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    scal(n, beta, y, 1);

    if (incX<0) {
        x -= incX*(m-1);
    }

    IntrinsicType  y_;
    IntrinsicType  A0_, A1_, A2_, A3_;
    IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
    IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

    IndexType j=0, jX=0;

    for (; j+4<m; j+=4, jX+=4*incX) {

        real_x0_.fill(real(alpha*x[jX       ]));
        imag_x0_.fill(imag(alpha*x[jX       ]));
        real_x1_.fill(real(alpha*x[jX+  incX]));
        imag_x1_.fill(imag(alpha*x[jX+  incX]));
        real_x2_.fill(real(alpha*x[jX+2*incX]));
        imag_x2_.fill(imag(alpha*x[jX+2*incX]));
        real_x3_.fill(real(alpha*x[jX+3*incX]));
        imag_x3_.fill(imag(alpha*x[jX+3*incX]));


        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );
            A1_.loadu(A+i+j*ldA+  ldA);
            A2_.loadu(A+i+j*ldA+2*ldA);
            A3_.loadu(A+i+j*ldA+3*ldA);


            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, real_x0_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A1_, real_x1_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A2_, real_x2_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A3_, real_x3_));


            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);


            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A0_, imag_x0_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A1_, imag_x1_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A2_, imag_x2_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A3_, imag_x3_));

            y_.storeu(y+i);
        }
    }


    for (; j<m; ++j, jX+=incX) {

        real_x0_.fill(real(alpha*x[jX  ]));
        imag_x0_.fill(imag(alpha*x[jX  ]));

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );

            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A0_, imag_x0_));

            y_.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);

    for (; i<n; ++i) {
        T y_;
        dotu(m, A+i, ldA, x, incX, y_);
        y[i] += alpha*y_;
    }

}


template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
gemv_complex_ct(IndexType m, IndexType n,
                const T &alpha,
                const T *A, IndexType ldA,
                const T *x, IndexType incX,
                const T &beta,
                T *y, IndexType incY)

{
    CXXBLAS_DEBUG_OUT("gemv_complex_ct");

    using std::imag;
    using std::real;

    ASSERT( incY==1 );

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    const int numElements = IntrinsicType::numElements;

    scal(n, beta, y, 1);

    if (incX<0) {
        x -= incX*(m-1);
    }

    IntrinsicType  y_;
    IntrinsicType  A0_, A1_, A2_, A3_;
    IntrinsicPrimitiveType real_x0_, real_x1_, real_x2_, real_x3_;
    IntrinsicPrimitiveType imag_x0_, imag_x1_, imag_x2_, imag_x3_;

    IndexType j=0, jX=0;

    for (; j+4<m; j+=4, jX+=4*incX) {

        real_x0_.fill(-real(alpha*x[jX       ]));
        imag_x0_.fill( imag(alpha*x[jX       ]));
        real_x1_.fill(-real(alpha*x[jX+  incX]));
        imag_x1_.fill( imag(alpha*x[jX+  incX]));
        real_x2_.fill(-real(alpha*x[jX+2*incX]));
        imag_x2_.fill( imag(alpha*x[jX+2*incX]));
        real_x3_.fill(-real(alpha*x[jX+3*incX]));
        imag_x3_.fill( imag(alpha*x[jX+3*incX]));


        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );
            A1_.loadu(A+i+j*ldA+  ldA);
            A2_.loadu(A+i+j*ldA+2*ldA);
            A3_.loadu(A+i+j*ldA+3*ldA);


            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A0_, real_x0_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A1_, real_x1_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A2_, real_x2_));
            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A3_, real_x3_));


            A0_ = intrinsic_swap_real_imag_(A0_);
            A1_ = intrinsic_swap_real_imag_(A1_);
            A2_ = intrinsic_swap_real_imag_(A2_);
            A3_ = intrinsic_swap_real_imag_(A3_);


            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, imag_x0_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A1_, imag_x1_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A2_, imag_x2_));
            y_ = intrinsic_add_(y_, intrinsic_mul_(A3_, imag_x3_));

            y_.storeu(y+i);
        }
    }


    for (; j<m; ++j, jX+=incX) {

        real_x0_.fill(-real(alpha*x[jX  ]));
        imag_x0_.fill( imag(alpha*x[jX  ]));

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            y_.loadu(y+i);

            A0_.loadu(A+i+j*ldA      );

            y_ = intrinsic_addsub_(y_, intrinsic_mul_(A0_, real_x0_));

            A0_ = intrinsic_swap_real_imag_(A0_);

            y_ = intrinsic_add_(y_, intrinsic_mul_(A0_, imag_x0_));

            y_.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);

    for (; i<n; ++i) {
        T y_;
        dot(m, A+i, ldA, x, incX, y_);
        y[i] += alpha*y_;
    }

}


template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value &&
                               flens::IsIntrinsicsCompatible<T>::value,
                               void>::Type
    gemv(StorageOrder order, Transpose transA,
         IndexType m, IndexType n,
         const T &alpha,
         const T *A, IndexType ldA,
         const T *x, IndexType incX,
         const T &beta,
         T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("gemv_intrinsics [real, " INTRINSIC_NAME "]");

    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        gemv(RowMajor, transA, n, m, alpha, A, ldA,
             x, incX, beta, y, incY);
        return;
    }

    if  ((transA==NoTrans || transA==Conj) && incX==1 ) {

        gemv_real_n(m, n, alpha, A, ldA, x, 1, beta, y, incY);

    } else if ((transA==Trans || transA==ConjTrans) && incY==1 ) {

        gemv_real_t(m, n, alpha, A, ldA, x, incX, beta, y, 1);

    } else {
        cxxblas::gemv<IndexType, T, T, T, T, T>(RowMajor, transA, m, n, alpha,
                                                A, ldA, x, incX, beta, y, incY);
    }
}

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsComplex<T>::value &&
                               flens::IsIntrinsicsCompatible<T>::value,
                               void>::Type
    gemv(StorageOrder order, Transpose transA,
         IndexType m, IndexType n,
         const T &alpha,
         const T *A, IndexType ldA,
         const T *x, IndexType incX,
         const T &beta,
         T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("gemv_intrinsics [complex, " INTRINSIC_NAME "]");

    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        gemv(RowMajor, transA, n, m, alpha, A, ldA,
             x, incX, beta, y, incY);
        return;
    }

    if  ( transA==NoTrans && incX==1 ) {

        gemv_complex_n(m, n, alpha, A, ldA, x, 1, beta, y, incY);

    } else if ( transA==Conj && incX==1 ) {

        gemv_complex_c(m, n, alpha, A, ldA, x, 1, beta, y, incY);

    } else if ( transA==Trans && incY==1 ) {

        gemv_complex_t(m, n, alpha, A, ldA, x, incX, beta, y, 1);

    } else if ( transA==ConjTrans && incY==1 ) {

        gemv_complex_ct(m, n, alpha, A, ldA, x, incX, beta, y, 1);

    } else {
        cxxblas::gemv<IndexType, T, T, T, T, T>(RowMajor, transA, m, n, alpha,
                                                A, ldA, x, incX, beta, y, incY);
    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_GEMV_TCC
