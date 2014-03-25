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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_TCC 1

#include <string.h>
#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T>
void
PackMatrixColToColMajor( IndexType k, const T *A, IndexType ldA, T *A_to, IndexType ldaA_to )
{

    switch(ldaA_to){
        case 2:      
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = a_ij_pntr[0];
                A_to[1] = a_ij_pntr[1];
                
                A_to += 2;
            }
            break;
        case 4:  
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = a_ij_pntr[0];
                A_to[1] = a_ij_pntr[1];
                A_to[2] = a_ij_pntr[2];
                A_to[3] = a_ij_pntr[3];
                
                A_to += 4;
            }
            break;
        case 8:  
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = a_ij_pntr[0];
                A_to[1] = a_ij_pntr[1];
                A_to[2] = a_ij_pntr[2];
                A_to[3] = a_ij_pntr[3];
                A_to[4] = a_ij_pntr[4];
                A_to[5] = a_ij_pntr[5];
                A_to[6] = a_ij_pntr[6];
                A_to[7] = a_ij_pntr[7];
                
                A_to += 8;
            }
            break;
        case 16:
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = a_ij_pntr[0];
                A_to[1] = a_ij_pntr[1];
                A_to[2] = a_ij_pntr[2];
                A_to[3] = a_ij_pntr[3];
                A_to[4] = a_ij_pntr[4];
                A_to[5] = a_ij_pntr[5];
                A_to[6] = a_ij_pntr[6];
                A_to[7] = a_ij_pntr[7];
                A_to[8] = a_ij_pntr[8];
                A_to[9] = a_ij_pntr[9];
                A_to[10] = a_ij_pntr[10];
                A_to[11] = a_ij_pntr[11];
                A_to[12] = a_ij_pntr[12];
                A_to[13] = a_ij_pntr[13];
                A_to[14] = a_ij_pntr[14];
                A_to[15] = a_ij_pntr[15];
                
                A_to += 16;
            }
            break;
        default: ASSERT(0);
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToColMajor_4( IndexType k, const T *A, IndexType ldA, T *A_to )
{

    for(IndexType j=0; j<k; j++){  /* loop over columns of A */
        const T *a_ij_pntr = A+j*ldA;
        
        *A_to++ = *(a_ij_pntr++);
        *A_to++ = *(a_ij_pntr++);
        *A_to++ = *(a_ij_pntr++);
        *A_to++ = *(a_ij_pntr++);
        
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToRowMajor( IndexType k, const T *B, IndexType ldB, T *b_to, IndexType ldB_to )
{
    const T * b_pntr[ldB_to];

    for (IndexType i=0; i<ldB_to; ++i) {
        b_pntr[i  ] = B+ i   *ldB;
    }

    for(IndexType i=0; i<k; ++i){  /* loop over rows of B */
        for (IndexType j=0; j<ldB_to; ++j) {
            *b_to++ = *b_pntr[j  ]++;
        }
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToRowMajor_4( IndexType k, const T *B, IndexType ldB, T *b_to )
{
    const T
    *b_i0_pntr = B,
    *b_i1_pntr = B+  ldB,
    *b_i2_pntr = B+2*ldB,
    *b_i3_pntr = B+3*ldB;

    for(IndexType i=0; i<k; i++){  /* loop over rows of B */
        *b_to++ = *b_i0_pntr++;
        *b_to++ = *b_i1_pntr++;
        *b_to++ = *b_i2_pntr++;
        *b_to++ = *b_i3_pntr++;
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToColMajor_conj( IndexType k, const T *A, IndexType ldA, T *A_to, IndexType ldA_to )
{

    switch(ldA_to){
        case 2:      
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = conjugate(a_ij_pntr[0]);
                A_to[1] = conjugate(a_ij_pntr[1]);
                
                A_to += 2;
            }
            break;
        case 4:  
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = conjugate(a_ij_pntr[0]);
                A_to[1] = conjugate(a_ij_pntr[1]);
                A_to[2] = conjugate(a_ij_pntr[2]);
                A_to[3] = conjugate(a_ij_pntr[3]);
                
                A_to += 4;
            }
            break;
        case 8:  
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = conjugate(a_ij_pntr[0]);
                A_to[1] = conjugate(a_ij_pntr[1]);
                A_to[2] = conjugate(a_ij_pntr[2]);
                A_to[3] = conjugate(a_ij_pntr[3]);
                A_to[4] = conjugate(a_ij_pntr[4]);
                A_to[5] = conjugate(a_ij_pntr[5]);
                A_to[6] = conjugate(a_ij_pntr[6]);
                A_to[7] = conjugate(a_ij_pntr[7]);
                
                A_to += 8;
            }
            break;
        case 16:
            for(IndexType j=0; j<k; j++){  /* loop over columns of A */
                const T *a_ij_pntr = A+j*ldA;
        
                A_to[0] = conjugate(a_ij_pntr[0]);
                A_to[1] = conjugate(a_ij_pntr[1]);
                A_to[2] = conjugate(a_ij_pntr[2]);
                A_to[3] = conjugate(a_ij_pntr[3]);
                A_to[4] = conjugate(a_ij_pntr[4]);
                A_to[5] = conjugate(a_ij_pntr[5]);
                A_to[6] = conjugate(a_ij_pntr[6]);
                A_to[7] = conjugate(a_ij_pntr[7]);
                A_to[8] = conjugate(a_ij_pntr[8]);
                A_to[9] = conjugate(a_ij_pntr[9]);
                A_to[10] = conjugate(a_ij_pntr[10]);
                A_to[11] = conjugate(a_ij_pntr[11]);
                A_to[12] = conjugate(a_ij_pntr[12]);
                A_to[13] = conjugate(a_ij_pntr[13]);
                A_to[14] = conjugate(a_ij_pntr[14]);
                A_to[15] = conjugate(a_ij_pntr[15]);
                
                A_to += 16;
            }
            break;
        default: ASSERT(0);
    }
    
    for(IndexType j=0; j<k; j++){  
        const T *a_ij_pntr = A+j*ldA;
        for (IndexType l=0; l<ldA_to; ++l) {
        	*(A_to+l  ) = conjugate(*(a_ij_pntr+l  ));
        }

        A_to += ldA_to;
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToColMajor_4_conj( IndexType k, const T *A, IndexType ldA, T *A_to )
{

    for(IndexType j=0; j<k; j++){  /* loop over columns of A */
        const T *a_ij_pntr = A+j*ldA;
        
        *A_to++ = conjugate(*(a_ij_pntr++));
        *A_to++ = conjugate(*(a_ij_pntr++));
        *A_to++ = conjugate(*(a_ij_pntr++));
        *A_to++ = conjugate(*(a_ij_pntr++));
        
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToRowMajor_conj( IndexType k, const T *B, IndexType ldB, T *b_to, IndexType ldB_to )
{

    const T * b_pntr[ldB_to];

    for (IndexType i=0; i<ldB_to; ++i) {
        b_pntr[i  ] = B+ i   *ldB;
    }

    for(IndexType i=0; i<k; ++i){  
        for (IndexType j=0; j<ldB_to; ++j) {
            *b_to++ = conjugate(*b_pntr[j  ]++);
        }
    }
}

template <typename IndexType, typename T>
void
PackMatrixColToRowMajor_4_conj( IndexType k, const T *B, IndexType ldB, T *b_to )
{
    const T
    *b_i0_pntr = B,
    *b_i1_pntr = B+  ldB,
    *b_i2_pntr = B+2*ldB,
    *b_i3_pntr = B+3*ldB;

    for(IndexType i=0; i<k; i++){  /* loop over rows of B */
        *b_to++ = conjugate(*b_i0_pntr++);
        *b_to++ = conjugate(*b_i1_pntr++);
        *b_to++ = conjugate(*b_i2_pntr++);
        *b_to++ = conjugate(*b_i3_pntr++);
    }
}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_PACKMATRIX_TCC
