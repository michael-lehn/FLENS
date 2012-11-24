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
    
        IntrinsicType _A0, _A1, _A2, _A3, _A4;
        IntrinsicType _x;
        IntrinsicType _y0, _y1, _y2, _y3;
        
        _y0.setZero();
        _y1.setZero();
        _y2.setZero();
        _y3.setZero();

        
        IndexType j=0;
        
        for (; j+numElements-1<n; j+=numElements) {
            _A0.loadu(A+i*ldA+j      );
            _A1.loadu(A+i*ldA+j+  ldA);
            _A2.loadu(A+i*ldA+j+2*ldA);
            _A3.loadu(A+i*ldA+j+3*ldA);
                                  
            _x.loadu(x+j);
            _y0 = _intrinsic_add(_y0, _intrinsic_mul(_A0, _x));
            _y1 = _intrinsic_add(_y1, _intrinsic_mul(_A1, _x));
            _y2 = _intrinsic_add(_y2, _intrinsic_mul(_A2, _x));
            _y3 = _intrinsic_add(_y3, _intrinsic_mul(_A3, _x));
            
        }
        
        _y0.storeu(&tmp_result[0][0]);
        _y1.storeu(&tmp_result[1][0]);
        _y2.storeu(&tmp_result[2][0]);
        _y3.storeu(&tmp_result[3][0]);
        
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
    
        IntrinsicType _A, _x, _y;
        _y.setZero();
        
        IndexType j=0;
        
        for (; j+numElements-1<n; j+=numElements) {
            _A.loadu(A+i*ldA+j);
            _x.loadu(x+j);
            _y = _intrinsic_add(_y, _intrinsic_mul(_A, _x));
        }
        T tmp_result[numElements];
        _y.storeu(tmp_result);
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
    
    
    IntrinsicType _y;
    IntrinsicType _beta(beta);
            
    IntrinsicType _A0, _A1, _A2, _A3;  
            
    IntrinsicType _x0, _x1, _x2, _x3;
          
    scal(n, beta, y, 1);
    
    if (incX<0) {
        x -= incX*(m-1);
    }
    
    IndexType j=0, jX=0;

    for (; j+3<m; j+=4, jX+=4*incX) {
            
        _x0.fill(alpha*x[jX       ]);
        _x1.fill(alpha*x[jX+1*incX]);
        _x2.fill(alpha*x[jX+2*incX]);
        _x3.fill(alpha*x[jX+3*incX]);

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );
            _A1.loadu(A+i+j*ldA+  ldA);
            _A2.loadu(A+i+j*ldA+2*ldA);
            _A3.loadu(A+i+j*ldA+3*ldA);
                    
            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _x0));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A1, _x1));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A2, _x2));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A3, _x3));

            _y.storeu(y+i);
        }
    }
            
    for (; j<m; ++j, jX+=incX) {
            
        _x0.fill(alpha*x[jX]);                   

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );

            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _x0));
                    
            _y.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);
    
    for (; i<n; ++i) {
        T _y;
        dot(m, A+i, ldA, x, incX, _y);
        y[i] += alpha*_y;
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
    
    IntrinsicType _A0, _A1, _A2, _A3;
    IntrinsicType  _x;
    IntrinsicType  _y0, _y1, _y2, _y3;
    IntrinsicPrimitiveType _real_x, _imag_x;
                    
    T tmp_result[8][numElements];

    IndexType i=0, iY=0;
    for (; i+3<m; i+=4, iY+=4*incY) {
                

    _y0.setZero();
    _y1.setZero();
    _y2.setZero();
    _y3.setZero();
                    
    IndexType j=0;
    for (; j+numElements-1<n; j+=numElements) {

        _A0.loadu(A+i*ldA+j);
        _A1.loadu(A+i*ldA+j+ldA);
        _A2.loadu(A+i*ldA+j+2*ldA);
        _A3.loadu(A+i*ldA+j+3*ldA);
                        
        _x.loadu(x+j);

        _real_x = _intrinsic_real(_x);
        _imag_x = _intrinsic_imag(_x);
            
        _y0 = _intrinsic_add(_y0, _intrinsic_mul(_A0, _real_x));
        _y1 = _intrinsic_add(_y1, _intrinsic_mul(_A1, _real_x));
        _y2 = _intrinsic_add(_y2, _intrinsic_mul(_A2, _real_x));
        _y3 = _intrinsic_add(_y3, _intrinsic_mul(_A3, _real_x));
                                                
                      
        _A0 = _intrinsic_swap_real_imag(_A0);
        _A1 = _intrinsic_swap_real_imag(_A1);
        _A2 = _intrinsic_swap_real_imag(_A2);
        _A3 = _intrinsic_swap_real_imag(_A3);
                        
        _y0 = _intrinsic_addsub(_y0, _intrinsic_mul(_A0, _imag_x));
        _y1 = _intrinsic_addsub(_y1, _intrinsic_mul(_A1, _imag_x));
        _y2 = _intrinsic_addsub(_y2, _intrinsic_mul(_A2, _imag_x));
        _y3 = _intrinsic_addsub(_y3, _intrinsic_mul(_A3, _imag_x));
                        
        }

        _y0.storeu(&tmp_result[0][0]);
        _y1.storeu(&tmp_result[1][0]);
        _y2.storeu(&tmp_result[2][0]);
        _y3.storeu(&tmp_result[3][0]);
                    
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

        _y0.setZero();
        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            _A0.loadu(A+i*ldA+j);
            _x.loadu(x+j);

            _real_x = _intrinsic_real(_x);
            _imag_x = _intrinsic_imag(_x);
             
            
            _y0 = _intrinsic_add(_y0, _intrinsic_mul(_A0, _real_x));

            _A0 = _intrinsic_swap_real_imag(_A0);

            _y0 = _intrinsic_addsub(_y0, _intrinsic_mul(_A0, _imag_x));
        }
                    
        _y0.storeu(&tmp_result[0][0]);
                    
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
    
    IntrinsicType _A0, _A1, _A2, _A3;
    IntrinsicType  _x;
    IntrinsicType  _y0, _y1, _y2, _y3;
    IntrinsicPrimitiveType _real_x, _imag_x, _minusOne;
    
    PT minusOne(-1);
    _minusOne.fill(minusOne);
                    
    T tmp_result[4][numElements];

    IndexType i=0, iY=0;
    for (; i+3<m; i+=4, iY+=4*incY) {
                

        _y0.setZero();
        _y1.setZero();
        _y2.setZero();
        _y3.setZero();
                    
        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            _A0.loadu(A+i*ldA+j);
            _A1.loadu(A+i*ldA+j+ldA);
            _A2.loadu(A+i*ldA+j+2*ldA);
            _A3.loadu(A+i*ldA+j+3*ldA);
                        
            _x.loadu(x+j);
        
            _real_x = _intrinsic_mul(_minusOne, _intrinsic_real(_x));
            _imag_x = _intrinsic_imag(_x);
            
            _y0 = _intrinsic_addsub(_y0, _intrinsic_mul(_A0, _real_x));
            _y1 = _intrinsic_addsub(_y1, _intrinsic_mul(_A1, _real_x));
            _y2 = _intrinsic_addsub(_y2, _intrinsic_mul(_A2, _real_x));
            _y3 = _intrinsic_addsub(_y3, _intrinsic_mul(_A3, _real_x));
                                                
                        
            _A0 = _intrinsic_swap_real_imag(_A0);
            _A1 = _intrinsic_swap_real_imag(_A1);
            _A2 = _intrinsic_swap_real_imag(_A2);
            _A3 = _intrinsic_swap_real_imag(_A3);
                        
            _y0 = _intrinsic_add(_y0, _intrinsic_mul(_A0, _imag_x));
            _y1 = _intrinsic_add(_y1, _intrinsic_mul(_A1, _imag_x));
            _y2 = _intrinsic_add(_y2, _intrinsic_mul(_A2, _imag_x));
            _y3 = _intrinsic_add(_y3, _intrinsic_mul(_A3, _imag_x));
                        
        }

        _y0.storeu(&tmp_result[0][0]);
        _y1.storeu(&tmp_result[1][0]);
        _y2.storeu(&tmp_result[2][0]);
        _y3.storeu(&tmp_result[3][0]);
                    
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

        _y0.setZero();
        IndexType j=0;
        for (; j+numElements-1<n; j+=numElements) {

            _A0.loadu(A+i*ldA+j);
            _x.loadu(x+j);

            _real_x = _intrinsic_mul(_minusOne, _intrinsic_real(_x));
            _imag_x = _intrinsic_imag(_x);
             
            
            _y0 = _intrinsic_addsub(_y0, _intrinsic_mul(_A0, _real_x));

            _A0 = _intrinsic_swap_real_imag(_A0);

            _y0 = _intrinsic_add(_y0, _intrinsic_mul(_A0, _imag_x));
        }
                    
        _y0.storeu(&tmp_result[0][0]);
                    
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
    
    IntrinsicType  _y;
    IntrinsicType  _A0, _A1, _A2, _A3;
    IntrinsicPrimitiveType _real_x0, _real_x1, _real_x2, _real_x3;
    IntrinsicPrimitiveType _imag_x0, _imag_x1, _imag_x2, _imag_x3;

    IndexType j=0, jX=0;

    for (; j+4<m; j+=4, jX+=4*incX) {
            
        _real_x0.fill(real(alpha*x[jX       ]));
        _imag_x0.fill(imag(alpha*x[jX       ]));
        _real_x1.fill(real(alpha*x[jX+  incX]));
        _imag_x1.fill(imag(alpha*x[jX+  incX]));
        _real_x2.fill(real(alpha*x[jX+2*incX]));
        _imag_x2.fill(imag(alpha*x[jX+2*incX]));
        _real_x3.fill(real(alpha*x[jX+3*incX]));
        _imag_x3.fill(imag(alpha*x[jX+3*incX]));
                

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );
            _A1.loadu(A+i+j*ldA+  ldA);
            _A2.loadu(A+i+j*ldA+2*ldA);
            _A3.loadu(A+i+j*ldA+3*ldA);


            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _real_x0));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A1, _real_x1));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A2, _real_x2));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A3, _real_x3));

                    
            _A0 = _intrinsic_swap_real_imag(_A0);
            _A1 = _intrinsic_swap_real_imag(_A1);
            _A2 = _intrinsic_swap_real_imag(_A2);
            _A3 = _intrinsic_swap_real_imag(_A3);


            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A0, _imag_x0));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A1, _imag_x1));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A2, _imag_x2));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A3, _imag_x3));
                        
            _y.storeu(y+i);
        }
    }
    

    for (; j<m; ++j, jX+=incX) {
            
        _real_x0.fill(real(alpha*x[jX  ]));
        _imag_x0.fill(imag(alpha*x[jX  ]));
                    
        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );

            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _real_x0));
                    
            _A0 = _intrinsic_swap_real_imag(_A0);
                        
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A0, _imag_x0));
                        
            _y.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);

    for (; i<n; ++i) {
        T _y;
        dotu(m, A+i, ldA, x, incX, _y);
        y[i] += alpha*_y;
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
    
    IntrinsicType  _y;
    IntrinsicType  _A0, _A1, _A2, _A3;
    IntrinsicPrimitiveType _real_x0, _real_x1, _real_x2, _real_x3;
    IntrinsicPrimitiveType _imag_x0, _imag_x1, _imag_x2, _imag_x3;

    IndexType j=0, jX=0;

    for (; j+4<m; j+=4, jX+=4*incX) {
            
        _real_x0.fill(-real(alpha*x[jX       ]));
        _imag_x0.fill( imag(alpha*x[jX       ]));
        _real_x1.fill(-real(alpha*x[jX+  incX]));
        _imag_x1.fill( imag(alpha*x[jX+  incX]));
        _real_x2.fill(-real(alpha*x[jX+2*incX]));
        _imag_x2.fill( imag(alpha*x[jX+2*incX]));
        _real_x3.fill(-real(alpha*x[jX+3*incX]));
        _imag_x3.fill( imag(alpha*x[jX+3*incX]));
                

        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );
            _A1.loadu(A+i+j*ldA+  ldA);
            _A2.loadu(A+i+j*ldA+2*ldA);
            _A3.loadu(A+i+j*ldA+3*ldA);


            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A0, _real_x0));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A1, _real_x1));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A2, _real_x2));
            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A3, _real_x3));

                    
            _A0 = _intrinsic_swap_real_imag(_A0);
            _A1 = _intrinsic_swap_real_imag(_A1);
            _A2 = _intrinsic_swap_real_imag(_A2);
            _A3 = _intrinsic_swap_real_imag(_A3);


            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _imag_x0));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A1, _imag_x1));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A2, _imag_x2));
            _y = _intrinsic_add(_y, _intrinsic_mul(_A3, _imag_x3));
                        
            _y.storeu(y+i);
        }
    }
    

    for (; j<m; ++j, jX+=incX) {
            
        _real_x0.fill(-real(alpha*x[jX  ]));
        _imag_x0.fill( imag(alpha*x[jX  ]));
                    
        for (IndexType i=0, ii=0; i+numElements-1<n; i+=numElements, ++ii) {
            _y.loadu(y+i);
                    
            _A0.loadu(A+i+j*ldA      );

            _y = _intrinsic_addsub(_y, _intrinsic_mul(_A0, _real_x0));
                    
            _A0 = _intrinsic_swap_real_imag(_A0);
                        
            _y = _intrinsic_add(_y, _intrinsic_mul(_A0, _imag_x0));
                        
            _y.storeu(y+i);
        }
    }

    IndexType i=n-(n%numElements);

    for (; i<n; ++i) {
        T _y;
        dot(m, A+i, ldA, x, incX, _y);
        y[i] += alpha*_y;
    }

}


template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsReal<T>::value, void>::Type
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
        cxxblas::gemv<IndexType, T, T, T, T, T>(RowMajor, transA, m, n, alpha, A, ldA, x, incX, beta, y, incY);
    }    
}

template <typename IndexType, typename T>
    typename flens::RestrictTo<flens::IsComplex<T>::value, void>::Type
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
        cxxblas::gemv<IndexType, T, T, T, T, T>(RowMajor, transA, m, n, alpha, A, ldA, x, incX, beta, y, incY);
    }    
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL2_GEMV_TCC
