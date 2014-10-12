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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_RACXPY_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_RACXPY_TCC 1

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
racxpy(IndexType n, const T &alpha, const T *x,
      IndexType incX, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("racxpy_intrinsics [real, " INTRINSIC_NAME "]");

    raxpy(n, alpha, x, incX, y, incY);
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value &&
                           flens::IsIntrinsicsCompatible<T>::value,
                           void>::Type
racxpy(IndexType n, const T &alpha, const T *x,
      IndexType incX, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("racxpy_intrinsics [complex, " INTRINSIC_NAME "]");

    using std::real;
    using std::imag;

    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;

    if (incX==1 && incY==1) {

        const int numElements = IntrinsicType::numElements;

        IndexType i=0;

        IntrinsicType x_, y_, tmp_;

        if (abs(real(alpha)) < abs(imag(alpha))) {

            PT r   = real(alpha)/imag(alpha);
            PT den = imag(alpha) + r*real(alpha);

            IntrinsicPrimitiveType r_(r);
            IntrinsicPrimitiveType den_(den);

           for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                y_.loadu(y+i);

                tmp_ = intrinsic_swap_real_imag_(x_);
                tmp_ = intrinsic_addsub_(tmp_, intrinsic_mul_(r_,x_));
                y_   = intrinsic_sub_(y_, intrinsic_div_(tmp_, den_));

                y_.storeu(y+i);
            }

        } else {

            PT r   = imag(alpha)/real(alpha);
            PT den = real(alpha) + r*imag(alpha);

            IntrinsicPrimitiveType r_(r);
            IntrinsicPrimitiveType den_(den);

           for (; i+numElements-1<n; i+=numElements) {
                x_.loadu(x+i);
                y_.loadu(y+i);

                tmp_ = intrinsic_mul_(r_,intrinsic_swap_real_imag_(x_));
                x_   = intrinsic_addsub_(tmp_, x_);
                y_   = intrinsic_sub_(y_, intrinsic_div_(x_, den_));

                y_.storeu(y+i);
            }



        }

        for (; i<n; ++i) {
            y[i] += conjugate(x[i])/alpha;
        }

    } else {

        cxxblas::racxpy<IndexType, T, T ,T>(n, alpha, x, incX, y, incY);

    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_RACXPY_TCC
