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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_AXPBY_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_AXPBY_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/level1/axpy.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
axpby(IndexType n, const T &alpha, const T *x,
      IndexType incX, const T &beta, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("axpby_intrinsics [real, " INTRINSIC_NAME "]");

//
//  Catch simple cases
//  1. y = beta*y
//  2. y = y + alpha*x
//  3. y = x
//
    if (alpha==T(0)) {
        cxxblas::scal(n, beta, y, incY);
        return;
    }
    if (beta==T(1)) {
        cxxblas::axpy(n, alpha, x, incX, y, incY);
    return;
    }
    if (alpha==T(1) && beta==T(0)) {
        cxxblas::copy(n, x, incX, y, incY);
    return;
    }

    if (incX==1 && incY==1) {
        typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
        const int numElements = IntrinsicType::numElements;

        IndexType i=0;

        IntrinsicType x_, y_;
        IntrinsicType alpha_(alpha);
        IntrinsicType beta_(beta);

//
//      Case y = alpha*x
//
        if (beta==T(0)) {
            for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                x_ = intrinsic_mul_(alpha_, x_);
                x_.storeu(y+i);
            }
        } else {
//
//          Case y = beta*y + alpha*x
//
            for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                y_.loadu(y+i);
                y_ = intrinsic_add_(intrinsic_mul_(beta_,y_),
                                    intrinsic_mul_(alpha_, x_));
                y_.storeu(y+i);
            }
        }
        for (; i<n; ++i) {
            y[i] = beta*y[i] + alpha*x[i];
        }

    } else {

        cxxblas::axpby<IndexType, T, T, T ,T>(n, alpha, x, incX, beta, y, incY);

    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
axpby(IndexType n, const T &alpha, const T *x,
     IndexType incX, const T &beta, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("axpby_intrinsics [complex, " INTRINSIC_NAME "]");

    using std::real;
    using std::imag;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

//
//  Catch simple cases
//  1. y = beta*y
//  2. y = y + alpha*x
//  3. y = x
//
    if (alpha==T(0)) {
        cxxblas::scal(n, beta, y, incY);
        return;
    }
    if (beta==T(1)) {
        cxxblas::axpy(n, alpha, x, incX, y, incY);
    return;
    }
    if (alpha==T(1) && beta==T(0)) {
        cxxblas::copy(n, x, incX, y, incY);
    return;
    }

//
//  Case y = beta*y + alpha*x
//
    if (incX==1 && incY==1) {


        if (imag(alpha)==PT(0) && imag(beta)==PT(0)) {
//
//          Only real scalars
//
            axpby(2*n, real(alpha), reinterpret_cast<const PT*>(x), 1,
                       real(beta) , reinterpret_cast<PT*>(y), 1);
           return;
        }

        IndexType i=0;

        if (beta==T(0)) {
//
//          beta is zero
//
            const int numElements = IntrinsicType::numElements;

            IntrinsicType x_, y_, tmp_;
            IntrinsicPrimitiveType real_alpha_(real(alpha));
            IntrinsicPrimitiveType imag_alpha_(imag(alpha));


            for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                y_ = intrinsic_mul_(real_alpha_, x_);
                x_ = intrinsic_swap_real_imag_(x_);
                y_ = intrinsic_addsub_(y_, intrinsic_mul_(imag_alpha_, x_));
                y_.storeu(y+i);
            }


        } else {
//
//          alpha and beta are complex
//
            const int numElements = IntrinsicType::numElements;

            IntrinsicType x_, y_, tmp_;
            IntrinsicPrimitiveType real_alpha_(real(alpha));
            IntrinsicPrimitiveType imag_alpha_(imag(alpha));
            IntrinsicPrimitiveType real_beta_(real(beta));
            IntrinsicPrimitiveType imag_beta_(imag(beta));

            for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                y_.loadu(y+i);

                tmp_ = intrinsic_mul_(real_beta_, y_);
                y_ = intrinsic_swap_real_imag_(y_);
                y_ = intrinsic_mul_(imag_beta_, y_);
                y_ = intrinsic_addsub_(tmp_, y_);

                y_ = intrinsic_add_(y_, intrinsic_mul_(real_alpha_, x_));
                x_ = intrinsic_swap_real_imag_(x_);
                y_ = intrinsic_addsub_(y_, intrinsic_mul_(imag_alpha_, x_));
                y_.storeu(y+i);
            }
        }
        for (; i<n; ++i) {
            y[i] = beta*y[i] + alpha*x[i];
        }

    } else {

        cxxblas::axpby<IndexType, T, T, T ,T>(n, alpha, x, incX, beta, y, incY);

    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_AXPBY_TCC
