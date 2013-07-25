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

#ifndef PLAYGROUND_CXXBLAS_LEVEL3EXTENSION_GEMM_STASSEN_TCC
#define PLAYGROUND_CXXBLAS_LEVEL3EXTENSION_GEMM_STASSEN_TCC 1

#include <cxxblas/cxxblas.h>

namespace cxxblas {


#ifdef USE_STRASSEN

#ifdef USE_COMPLEX_3M

template <typename IndexType, typename T>
void
gemm_complex_3m_kernel(StorageOrder order,
                       Transpose transA, Transpose transB,
                       IndexType m, IndexType n, IndexType k,
                       T alpha,
                       const T *A, IndexType ldA,
                       const T *B, IndexType ldB,
                       T beta,
                       T *C, IndexType ldC,
                       T *ptmp = NULL,
                       size_t size_tmp = 0)
{
/*
 *  Function for complex matrix-matrix multiplication
 *  with only three real-valued multiplications.
 *  Is called after the Strassen-Algorithm
 */

    typedef typename flens::ComplexTrait<T>::PrimitiveType PT;

    PT one(1), zero(0);
    ASSERT(2*size_tmp>=size_t(m*k)+size_t(k*n));
  
    // Temporary objects for real/imaginary parts of A and B
    PT *A_tmp = reinterpret_cast<PT*>(ptmp);
    PT *B_tmp = A_tmp + m*k;

    IndexType ld_Atmp = (transA==NoTrans || transA==Conj) ? m : k;
    IndexType ld_Btmp = (transB==NoTrans || transB==Conj) ? k : n;
 
    PT mod_one_A = (transA==NoTrans || transA==Trans) ? PT(1) : PT(-1);
    PT mod_one_B = (transB==NoTrans || transB==Trans) ? PT(1) : PT(-1);

    Transpose mod_transA = (transA==NoTrans || transA==Conj) ? NoTrans : Trans;
    Transpose mod_transB = (transB==NoTrans || transB==Conj) ? NoTrans : Trans;

    PT r_alpha = (cxxblas::imag(alpha==PT(0))) ? cxxblas::real(alpha) : PT(1);
    ASSERT(beta==T(0));

    // A_tmp = A_r + A_i
    if (transA==NoTrans || transA==Conj) {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(m,
                          reinterpret_cast<const PT*>(A+i*ldA),IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
            cxxblas::axpy(m, mod_one_A, 
                          reinterpret_cast<const PT*>(A+i*ldA)+1,IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            cxxblas::copy(k, 
                          reinterpret_cast<const PT*>(A+i*ldA), IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
            cxxblas::axpy(k, mod_one_A, 
                          reinterpret_cast<const PT*>(A+i*ldA)+1,IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    }
    // B_tmp = B_r
    if (transB==NoTrans || transB==Conj) {
        for (IndexType i=0; i<n; ++i) {
            cxxblas::copy(k,
                          reinterpret_cast<const PT*>(B+i*ldB), IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(n, 
                          reinterpret_cast<const PT*>(B+i*ldB),IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    }

    // C_up = r_alpha*A_tmp*B_tmp
    cxxblas::gemm(order, mod_transA, mod_transB, 
                  m, n, k, 
                  r_alpha, 
                  A_tmp, ld_Atmp, 
                  B_tmp, ld_Btmp, 
                  zero, 
                  reinterpret_cast<PT *>(C), 2*ldC);


    //A_tmp = A_r
    if (transA==NoTrans) {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(m, 
                          reinterpret_cast<const PT*>(A+i*ldA), IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            cxxblas::copy(k, 
                          reinterpret_cast<const PT*>(A+i*ldA), IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    }
    // B_tmp = B_r-B_i
    if (transB==NoTrans) {
        for (IndexType i=0; i<n; ++i) {
            cxxblas::copy(k, 
                          reinterpret_cast<const PT*>(B+i*ldB), IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
            cxxblas::axpy(k, -mod_one_B, 
                          reinterpret_cast<const PT*>(B+i*ldB)+1, IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(n, 
                          reinterpret_cast<const PT*>(B+i*ldB), IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
            cxxblas::axpy(n, -mod_one_B, 
                          reinterpret_cast<const PT*>(B+i*ldB)+1, IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    }

    // C_down = r_alpha*A_tmp*B_tmp
    cxxblas::gemm(order, mod_transA, mod_transB, 
                  m, n, k, 
                  r_alpha, 
                  A_tmp, ld_Atmp, 
                  B_tmp, ld_Btmp, 
                  zero, 
                  reinterpret_cast<PT *>(C)+m, 2*ldC);

    // C_up = C_up - C_down
    cxxblas::geaxpy(order, NoTrans, 
                    m, n, 
                   -one, 
                   reinterpret_cast<PT *>(C)+m, 2*ldC, 
                   reinterpret_cast<PT *>(C), 2*ldC);

    // A_tmp = A_r - A_i
    if (transA==NoTrans || transA==Conj) {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(m, 
                          reinterpret_cast<const PT*>(A+i*ldA), IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
            cxxblas::axpy(m, -mod_one_A, 
                          reinterpret_cast<const PT*>(A+i*ldA)+1, IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            cxxblas::copy(k, 
                          reinterpret_cast<const PT*>(A+i*ldA), IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
            cxxblas::axpy(k, -mod_one_A, 
                          reinterpret_cast<const PT*>(A+i*ldA)+1, IndexType(2),
                          A_tmp+i*ld_Atmp, IndexType(1));
        }
    }

    // B_tmp = B_i
    if (transB==NoTrans) {
        for (IndexType i=0; i<n; ++i) {
            cxxblas::copy(k, 
                          reinterpret_cast<const PT*>(B+i*ldB)+1, IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    } else {
        for (IndexType i=0; i<k; ++i) {
            cxxblas::copy(n, 
                          reinterpret_cast<const PT*>(B+i*ldB)+1, IndexType(2),
                          B_tmp+i*ld_Btmp, IndexType(1));
        }
    }

    // C_down = C_down + A_tmp*B_tmp
    cxxblas::gemm(order, mod_transA, mod_transB, 
                  m, n, k, 
                  mod_one_B*r_alpha, 
                  A_tmp, ld_Atmp,
                  B_tmp, ld_Btmp, 
                  one, 
                  reinterpret_cast<PT *>(C)+m, 2*ldC);


    // Space for temporary Column of C
    PT *tmpC_row = A_tmp;
    if (k<2) {
        tmpC_row = new PT[2*m];
    }
    
    // Merge C = C_r + i*C_i
    for (IndexType i=0; i<n; ++i) {
        cxxblas::copy(m, 
                      C+i*ldC, IndexType(1), 
                      reinterpret_cast<std::complex<PT> *>(tmpC_row), IndexType(1));

	// copy real part
        cxxblas::copy(m, 
                      tmpC_row+m, IndexType(1), 
                      reinterpret_cast<PT *>(C+i*ldC), IndexType(2));
        // copy imaginary part
        cxxblas::copy(m, 
                      tmpC_row, IndexType(1), 
                      reinterpret_cast<PT *>(C+i*ldC)+1, IndexType(2));
        
       if (!cxxblas::imag(alpha==PT(0))) {
           cxxblas::scal(m, alpha, C+i*ldC, IndexType(1));
       }

    }
    

    if (k<2) {
        delete[] tmpC_row;
    }

}

#endif // USE_COMPLEX_3M

template <typename IndexType, typename T>
void
gemm_strassen_kernel(StorageOrder order,
                     Transpose transA, Transpose transB,
                     IndexType m, IndexType n, IndexType k,
                     T alpha,
                     const T *A, IndexType ldA,
                     const T *B, IndexType ldB,
                     T beta,
                     T *C, IndexType ldC,
                     T *ptmp = NULL,
                     size_t size_tmp = 0)
{
/*
 * Implementation of the Winograd Variant of Strassen's 
 * Matrix-Matrix Multiply Algorithm. It is based on
 * 
 * C Douglas et al: 
 * "GEMMW: A Portable Level 3 BLAS  Winograd Variant 
 *   of Strassen's Matrix-Matrix Multiply Algorithm."
 * J. Comput. Phys. 110, 1-10, 1994. 
 *
 */

    ASSERT(beta==T(0));

    T one(1), zero(0);
    IndexType min_mnk = std::min(m, std::min(n,k));

    if (min_mnk<=MINDIM) {
        ASSERT(size_tmp>=0);
#       ifndef USE_COMPLEX_3M
            cxxblas::gemm(order, transA, transB, 
                          m, n, k, 
                          alpha, 
                          A, ldA, 
                          B, ldB, 
                          beta, 
                          C, ldC);
	    return;
#       else 

            if(flens::IsReal<T>::value) {
                cxxblas::gemm(order, transA, transB, 
                              m, n, k, 
			      alpha, 
                              A, ldA, 
                              B, ldB, 
                              beta, 
                              C, ldC);
            } else {
                cxxblas::gemm_complex_3m_kernel(order, transA, transB, 
                                                m, n, k, 
                                                alpha, 
                                                A, ldA, 
                                                B, ldB, 
                                                beta, 
                                                C, ldC, 
                                                ptmp, size_tmp);
            }
            return;

#       endif

    }

    bool k_odd = !(k%2==0);

    IndexType mh =  m   /2;
    IndexType nh =  n   /2;
    IndexType kh = (k+1)/2;

    const T *pA11 = A;
    const T *pA21 = (transA==NoTrans || transA==Conj) ? A+mh        : A   +mh*ldA;
    const T *pA12 = (transA==NoTrans || transA==Conj) ? A   +kh*ldA : A+kh       ;
    const T *pA22 = (transA==NoTrans || transA==Conj) ? A+mh+kh*ldA : A+kh+mh*ldA;

    const T *pB11 = B;
    const T *pB21 = (transB==NoTrans || transB==Conj) ? B+kh        : B   +kh*ldB;
    const T *pB12 = (transB==NoTrans || transB==Conj) ? B   +nh*ldB : B+nh       ;
    const T *pB22 = (transB==NoTrans || transB==Conj) ? B+kh+nh*ldB : B+nh+kh*ldB;

    T *pC11 = C;
    T *pC21 = C + mh;
    T *pC12 = C      + nh*ldC;
    T *pC22 = C + mh + nh*ldC;

    IndexType A_numRows_h = (transA==NoTrans || transA==Conj) ? mh : kh;
    IndexType A_numCols_h = (transA==NoTrans || transA==Conj) ? kh : mh;

    IndexType B_numRows_h = (transB==NoTrans || transB==Conj) ? kh : nh;
    IndexType B_numCols_h = (transB==NoTrans || transB==Conj) ? nh : kh;


    // Temporary Storage
    T *ptmpB = ptmp;
    T *ptmpA = ptmp + nh*kh;
    T *ptmpC = ptmp + nh*kh;

    ptmp = ptmp + (nh*kh+mh*std::max(kh,nh));
    ASSERT(size_tmp >= size_t(nh*kh) + size_t(mh*std::max(kh,nh)));
    size_tmp = size_tmp - nh*kh - mh*std::max(kh,nh);

    // Leading Dimensions of temporary Storage
    IndexType ld_tmpA = A_numRows_h;
    IndexType ld_tmpB = B_numRows_h;
    IndexType ld_tmpC = mh;


    // Step 1: tmpB  = B22 - B12;
    if (!k_odd) {
        //
        cxxblas::gecopy(order, NoTrans, 
                        B_numRows_h, B_numCols_h, 
                        pB22, ldB, 
                        ptmpB, ld_tmpB);
        cxxblas::geaxpy(order, NoTrans, 
                        B_numRows_h, B_numCols_h, 
                        -one, 
                        pB12, ldB, 
                        ptmpB, ld_tmpB);
    } else {
        if (transB==NoTrans) {
            cxxblas::gecopy(order, NoTrans,
                            kh-1, nh, 
                            pB22, ldB, 
                            ptmpB, ld_tmpB);
            cxxblas::scal(nh, zero, ptmpB+kh-1, ld_tmpB);
            cxxblas::geaxpy(order, NoTrans, 
                            kh, nh, 
                            -one, 
                            pB12, ldB, 
                            ptmpB, ld_tmpB);
        } else {
            cxxblas::scal(nh, zero, ptmpB+(kh-1)*ld_tmpB, IndexType(1));
            cxxblas::gecopy(order, NoTrans, 
                            nh, kh-1, 
                            pB22, ldB, 
                            ptmpB, ld_tmpB);
            cxxblas::geaxpy(order, NoTrans, 
                            nh, kh,
                            -one,
                            pB12, ldB,
                            ptmpB, ld_tmpB);
        }
    }


    // Step 2: tmpA  = A11 - A21;
    cxxblas::gecopy(order, NoTrans, 
                    A_numRows_h, A_numCols_h, 
                    pA11, ldA, 
                    ptmpA, ld_tmpA);
 
    cxxblas::geaxpy(order, NoTrans, 
                    A_numRows_h, A_numCols_h,
                    -one,  
                    pA21,ldA, 
                    ptmpA, ld_tmpA);

    // Step 3: C21 = alpha*op(tmpA)*op(tmpB)
    cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                  mh, nh, kh, 
                                  alpha, 
                                  ptmpA, ld_tmpA, 
                                  ptmpB, ld_tmpB, 
                                  zero, 
                                  pC21, ldC, 
                                  ptmp, size_tmp);

    // Step 4: tmpA = A21 + A22;
    if (!k_odd) {
        cxxblas::gecopy(order, NoTrans, 
                        A_numRows_h, A_numCols_h,
                        pA21, ldA, 
                        ptmpA, ld_tmpA);
        cxxblas::geaxpy(order, NoTrans, 
                        A_numRows_h, A_numCols_h,
                        one, 
                        pA22, ldA, 
                        ptmpA, ld_tmpA);
    } else {
        if (transA==NoTrans || transA==Conj) {
            cxxblas::gecopy(order, NoTrans, 
                            A_numRows_h, A_numCols_h, 
                            pA21, ldA, 
                            ptmpA, ld_tmpA);
            cxxblas::geaxpy(order, NoTrans, 
                            A_numRows_h, A_numCols_h-1,
                            one,
                            pA22, ldA, 
                            ptmpA, ld_tmpA);
        } else {
            cxxblas::gecopy(order, NoTrans, 
                            A_numRows_h, A_numCols_h, 
                            pA21, ldA,
                            ptmpA, ld_tmpA);
            cxxblas::geaxpy(order, NoTrans, 
                            A_numRows_h-1, A_numCols_h,
                            one, 
                            pA22, ldA, 
                            ptmpA, ld_tmpA);
        }
    }

    // Step 5: tmpB = B12 - B11;
    cxxblas::gecopy(order, NoTrans, 
                    B_numRows_h, B_numCols_h,
                    pB12, ldB, 
                    ptmpB, ld_tmpB);
    cxxblas::geaxpy(order, NoTrans, 
                    B_numRows_h, B_numCols_h, 
                    -one, 
                    pB11, ldB, 
                    ptmpB, ld_tmpB);

    // Step 6: C22 = alpha*op(tmpA)*op(tmpB);
    cxxblas::gemm_strassen_kernel(order, transA, transB,
                                  mh, nh, kh, 
                                  alpha, 
                                  ptmpA, ld_tmpA, 
                                  ptmpB, ld_tmpB, 
                                  zero,  
                                  pC22, ldC, 
                                  ptmp, size_tmp);

    // Step 7: tmpB = tmpB - B22;
    if (!k_odd) {
        cxxblas::geaxpy(order, NoTrans, 
                        B_numRows_h, B_numCols_h, 
                        -one, 
                        pB22, ldB, 
                        ptmpB, ld_tmpB);
    } else {
        if (transB==NoTrans || transB==Conj) {
            cxxblas::geaxpy(order, NoTrans, 
                            B_numRows_h-1, B_numCols_h,
                            -one, 
                            pB22, ldB, 
                            ptmpB, ld_tmpB);
        } else {
            cxxblas::geaxpy(order, NoTrans,
                            B_numRows_h, B_numCols_h-1,
                            -one,
                            pB22, ldB, 
                            ptmpB, ld_tmpB);
        }
    }

    // Step 8: tmpA = tmpA - A11;
    cxxblas::geaxpy(order, NoTrans, 
                    A_numRows_h, A_numCols_h, 
                    -one, 
                    pA11, ldA, 
                    ptmpA, ld_tmpA);

    // Step 9: C11 = -alpha*op(tmpA)*op(tmpB);
    cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                  mh, nh, kh, 
                                  -alpha, 
                                  ptmpA, ld_tmpA, 
                                  ptmpB, ld_tmpB, 
                                  zero, 
                                  pC11, ldC, 
                                  ptmp, size_tmp);

    // Step 10: tmpA = tmpA - A12;
    if (!k_odd) {

        cxxblas::geaxpy(order, NoTrans, 
                        A_numRows_h, A_numCols_h, 
                        -one, 
                        pA12, ldA, 
                        ptmpA, ld_tmpA);

    } else {
        if (transA==NoTrans || transA==Conj) {

            cxxblas::geaxpy(order, NoTrans, 
                            A_numRows_h, A_numCols_h-1, 
                           -one, 
                           pA12, ldA, 
                           ptmpA, ld_tmpA);

        } else {

            cxxblas::geaxpy(order, NoTrans,
                            A_numRows_h-1, A_numCols_h, 
                           -one, 
                           pA12, ldA, 
                           ptmpA, ld_tmpA);
        }
    }

    // Step 11: C12 = -alpha*op(tmpA)*op(B22);
    if (!k_odd) {
        cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                      mh, nh, kh, 
                                      -alpha, 
                                      ptmpA, ld_tmpA, 
                                      pB22, ldB, 
                                      zero, 
                                      pC12, ldC, 
                                      ptmp, size_tmp);
    } else {
        if (transA==NoTrans || transA==Conj) {
                      
            cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                          mh, nh, kh-1,  
                                          -alpha, 
                                          ptmpA, ld_tmpA,
                                          pB22, ldB, 
                                          zero, 
                                          pC12, ldC, 
                                          ptmp, size_tmp);
        } else {
            cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                          mh, nh, kh-1, 
                                          -alpha, 
                                          ptmpA, ld_tmpA, 
                                          pB22, ldB, 
                                          zero, 
                                          pC12, ldC, 
                                          ptmp, size_tmp);
        }
    }

    // Step 12: C12 = C12 + C22;
    cxxblas::geaxpy(order, NoTrans, 
                    mh, nh, 
                    one, 
                    pC22, ldC,
                    pC12, ldC);

    // Step 13: tmpC = alpha*op(A11)*op(B11);
    cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                  mh, nh, kh, 
                                  alpha, 
                                  pA11, ldA, 
                                  pB11, ldB, 
                                  zero, 
                                  ptmpC, ld_tmpC, 
                                  ptmp, size_tmp);

    // Step 14: C11 = C11 + tmpC;
    cxxblas::geaxpy(order, NoTrans, 
                    mh, nh, 
                    one, 
                    ptmpC, ld_tmpC, 
                   pC11, ldC);

    // Step 15: C12 = C12 + C11;
    cxxblas::geaxpy(order, NoTrans, 
                    mh, nh,
                    one, 
                    pC11, ldC, 
                    pC12, ldC);

    // Step 16: C11 = C11 + C21;
    cxxblas::geaxpy(order, NoTrans, 
                    mh, nh, 
                    one, 
                    pC21, ldC, 
                    pC11, ldC);

    // Step 17: tmpB = tmpB + B21;
    if (!k_odd) {

        cxxblas::geaxpy(order, NoTrans, 
                        B_numRows_h, B_numCols_h, 
                        one, 
                        pB21, ldB,
                        ptmpB, ld_tmpB);

    } else {
        if (transB==NoTrans || transB==Conj) {

            cxxblas::geaxpy(order, NoTrans, 
                            B_numRows_h-1, B_numCols_h, 
                            one, 
                            pB21, ldB, 
                            ptmpB, ld_tmpB);

        } else {

            cxxblas::geaxpy(order, NoTrans, 
                            B_numRows_h, B_numCols_h-1,
                            one, 
                            pB21, ldB, 
                            ptmpB, ld_tmpB);
        }
    }

    // Step 18: C12 = -alpha*op(A22)*op(tmpB);
    if (!k_odd) {
        cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                      mh, nh, kh, 
                                      alpha, 
                                      pA22, ldA, 
                                      ptmpB, ld_tmpB, 
                                      zero, 
                                      pC21, ldC, 
                                      ptmp, size_tmp);
    } else {
        if (transB==NoTrans || transB==Conj) {
            cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                          mh, nh, kh-1, 
                                          alpha, 
                                          pA22, ldA, 
                                          ptmpB, ld_tmpB, 
                                          zero, 
                                          pC21, ldC, 
                                          ptmp, size_tmp);
        }  else {
            cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                          mh, nh, kh-1, 
                                          alpha, 
                                          pA22, ldA, 
                                          ptmpB, ld_tmpB, 
                                          zero, 
                                          pC21, ldC, 
                                          ptmp, size_tmp);
        }
    }
    // Step 19: C21 = C21  + C11;
    cxxblas::geaxpy(order, NoTrans,
                    mh, nh,
                    one, 
                    pC11, ldC, 
                    pC21, ldC);

    // Step 20: C22 = C22 + C11;
    cxxblas::geaxpy(order, NoTrans,   
                    mh, nh, 
                    one,
                    pC11, ldC, 
                    pC22, ldC);

    // Step 21: C11 = alpha*op(A12)*op(B21);
    if (!k_odd) {
        cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                      mh, nh, kh,
                                      alpha, 
                                      pA12, ldA, 
                                      pB21, ldB, 
                                      zero, 
                                      pC11, ldC, 
                                      ptmp, size_tmp);
    } else {
    	cxxblas::gemm_strassen_kernel(order, transA, transB,
                                      mh, nh, kh-1, 
                                      alpha, 
                                      pA12, ldA, 
                                      pB21, ldB, 
                                      zero, 
                                      pC11, ldC,
                                      ptmp, size_tmp);
    }

    // Step 22: C11 = C11 + tmpC;
    cxxblas::geaxpy(order, NoTrans, 
                    mh, nh, 
                    one,
                    ptmpC, ld_tmpC,
                    pC11, ldC);

}

template <typename IndexType, typename T>
void
gemm_strassen(StorageOrder order,
              Transpose transA, Transpose transB,
              IndexType m, IndexType n, IndexType k,
              T alpha,
              const T *A, IndexType ldA,
              const T *B, IndexType ldB,
              T beta,
              T *C, IndexType ldC)
{
        CXXBLAS_DEBUG_OUT("gemm strassen algorithm");

        using std::min;

        T one(1), zero(0);

        IndexType min_mnk = min(m,min(n,k));

        IndexType num_divides = 0;

        IndexType test_min_mnk = min_mnk;
        while (test_min_mnk>MINDIM) {
                test_min_mnk = (test_min_mnk+1)/2;
                ++num_divides;
        }

        IndexType mr = m - (m%IndexType(std::pow(IndexType(2),num_divides)));
        IndexType nr = n - (n%IndexType(std::pow(IndexType(2),num_divides)));

        num_divides = 0;

        test_min_mnk = min(k, min(mr, nr));
        while (test_min_mnk>MINDIM) {
                test_min_mnk = (test_min_mnk+1)/2;
                ++num_divides;
        }


        IndexType nh = nr/2;
        IndexType mh = mr/2;
        IndexType kh = (k+1)/2;

        // Calculate Storage
        size_t size_tmp = 0;
        for (IndexType i=1; i<=num_divides; ++i)
        {
                size_tmp += nh*kh+mh*std::max(kh,nh);
                mh = mh/2;
                nh = nh/2;
                kh = (kh+1)/2;
        }

#       ifdef USE_COMPLEX_3M
            // Additional Storage for three-sum trick
            if (flens::IsComplex<T>::value) {
                size_tmp += ((mh+1)*(kh+1))*2 +((kh+1)*(nh+1))*2;
            }
#       endif

        bool m_odd = !(m==mr);
        bool n_odd = !(n==nr);

        T *ptmp = new T[size_tmp];
        if (beta==zero) {
            cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                          mr, nr, k, 
                                          alpha, 
                                          A, ldA, 
                                          B, ldB, 
                                          zero, 
                                          C, ldC, 
                                          ptmp, size_tmp);
        } else {
                T *tmp_C = new T[mr*nr];
                cxxblas::gemm_strassen_kernel(order, transA, transB, 
                                              mr, nr, k, 
                                              alpha, 
                                              A, ldA, 
                                              B, ldB, 
                                              zero, 
                                              tmp_C, mr, 
                                              ptmp, size_tmp);
                cxxblas::gescal(order, mr, nr, beta, C, ldC);
                cxxblas::geaxpy(order, NoTrans, 
                                mr, nr,
                                one, 
                                tmp_C, mr, 
                                C, ldC);
                delete[] tmp_C;
        }
        delete[] ptmp;

        // Handle rows and columns with were stripped off
        if (n_odd && !m_odd) {
                if ( (transB==NoTrans) || (transB==Conj) ) {
                    cxxblas::gemm(order, transA, transB, 
                                  m, n-nr, k, 
                                  alpha, 
                                  A, ldA, 
                                  B+nr*ldB, ldB, 
                                  beta, 
                                  C+nr*ldC, ldC);
                } else {
                    cxxblas::gemm(order, transA, transB,
                                  m, n-nr, k,
                                  alpha, 
                                  A, ldA, 
                                  B+nr, ldB, 
                                  beta, 
                                  C+nr*ldC, ldC);
                }
        } else if (!n_odd && m_odd){
                if ( (transA==NoTrans) || (transA==Conj) ) {
                    cxxblas::gemm(order, transA, transB, 
                                  m-mr, n, k, 
                                  alpha, 
                                  A+mr, ldA, 
                                  B, ldB, 
                                  beta, 
                                  C+mr, ldC);
                } else {
                    cxxblas::gemm(order, transA, transB, 
                                  m-mr, n, k, 
                                  alpha, 
                                  A+mr*ldA, ldA, 
                                  B, ldB, 
                                  beta, 
                                  C+mr, ldC);
                }
        } else if (n_odd && m_odd) {
                if ( (transB==NoTrans) || (transB==Conj) ) {
                    cxxblas::gemm(order, transA, transB, 
                                  m, n-nr, k, 
                                  alpha, 
                                  A, ldA, 
                                  B+nr*ldB, ldB, 
                                  beta, 
                                  C+nr*ldC, ldC);
                } else {
                    cxxblas::gemm(order, transA, transB, 
                                  m, n-nr, k, 
                                  alpha,
                                  A, ldA, 
                                  B+nr, ldB, 
                                  beta, 
                                  C+nr*ldC, ldC);
                }
                if ( (transA==NoTrans) || (transA==Conj) ) {
                    cxxblas::gemm(order, transA, transB, 
                                  m-mr, nr, k,
                                  alpha, 
                                  A+mr, ldA, 
                                  B, ldB, 
                                  beta, 
                                  C+mr, ldC);
                } else {
                    cxxblas::gemm(order, transA, transB, 
                                  m-mr, nr, k, 
                                  alpha, 
                                  A+mr*ldA, ldA, 
                                  B, ldB, 
                                  beta, 
                                  C+mr, ldC);
                }
        }

        return;
}

#endif // USE_STRASSEN


} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_LEVEL3EXTENSION_GEMM_STASSEN_TCC
