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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_SCAL_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_SCAL_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
scal(IndexType n, const T &alpha, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("scal_intrinsics [real, " INTRINSIC_NAME "]");
    
    if (alpha==T(1))
        return;
        
    if (incY==1) {
        typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL> IntrinsicType;
        const int numElements = IntrinsicType::numElements;
        
        IndexType i=0;

        if (alpha==T(0)) {
                
            IntrinsicType _zero(T(0));
            for(;i+numElements-1<n;i+=numElements) {
                _zero.storeu(y+i);
            }
            
            for (;i<n;++i) {
                y[i] = T(0);
            }
            

        } else {
        

                
            IntrinsicType _alpha(alpha);
            IntrinsicType _y;
            for(;i+numElements-1<n;i+=numElements) {
                _y.loadu(y+i);
                _y = _intrinsic_mul(_alpha, _y);
                _y.storeu(y+i);
            }
            
            
            for (;i<n;++i) {
                y[i] *= alpha;
            }

        }
    } else {
        cxxblas::scal<IndexType, T, T>(n, alpha, y, incY);
    }
}

template <typename IndexType, typename T>
typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
scal(IndexType n, const T &alpha, T *y, IndexType incY)
{
    CXXBLAS_DEBUG_OUT("scal_intrinsics [complex, " INTRINSIC_NAME "]");

    using std::real;
    using std::imag;
    
    typedef Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>     IntrinsicType;
    typedef typename IntrinsicType::PrimitiveDataType  PT;
    typedef Intrinsics<PT, DEFAULT_INTRINSIC_LEVEL>    IntrinsicPrimitiveType;
    
    if (alpha==T(1))
        return;
        
    if (incY==1) {
    
        if (imag(alpha) == PT(0)) {
            scal(2*n, real(alpha), reinterpret_cast<PT*>(y), 1);
            return;
        }
        
        const int numElements = IntrinsicType::numElements;
        
        IndexType i=0;

        if (alpha==T(0)) {
                
            IntrinsicType _zero(PT(0));
            for(;i+numElements-1<n;i+=numElements) {
                _zero.storeu(y+i);
            }
            
            for (;i<n;++i) {
                y[i] = T(0);
            }

        } else {
                
            IntrinsicPrimitiveType _real_alpha(real(alpha));
            IntrinsicPrimitiveType _imag_alpha(imag(alpha));
            IntrinsicType _y, _tmp;
            
            for(;i+numElements-1<n;i+=numElements) {
                _y.loadu(y+i);
                _tmp = _intrinsic_mul(_real_alpha, _y);
                _y = _intrinsic_swap_real_imag(_y);
                _y = _intrinsic_mul(_imag_alpha, _y);
                _y = _intrinsic_addsub(_tmp, _y);
                _y.storeu(y+i);
            }
            
            for (;i<n;++i) {
                y[i] *= alpha;
            }

        }
    } else {
        cxxblas::scal<IndexType, T, T>(n, alpha, y, incY);
    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL1_SCAL_TCC
