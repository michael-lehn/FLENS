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

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_TCC 1

#include <cxxblas/cxxblas.h>
#include <playground/cxxblas/intrinsics/auxiliary/auxiliary.h>
#include <playground/cxxblas/intrinsics/includes.h>
#include <playground/cxxblas/intrinsics/level3/gemm/packmatrix.h>
#include <vector>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_n_n(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;

    for ( j=0; j+3<n; j+=4 ){
        if ( first_call )
            PackMatrixColToRowMajor_4( k, B+j*ldB, ldB, &B_packed[ j*k ] );

        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor(k, A+i, ldA, &A_packed[ i*k ], 2*numElements);
            }

            kernel_gemm_2numElementsx4(k, alpha, A_packed+i*k, 2*numElements, B_packed+j*k, k, beta, C+i+j*ldC, ldC);
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*B[l+ j   *ldB];
                tmp_j1 += A[i+l*ldA]*B[l+(j+1)*ldB];
                tmp_j2 += A[i+l*ldA]*B[l+(j+2)*ldB];
                tmp_j3 += A[i+l*ldA]*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i, ldA, B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }
}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_n_t(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*B[l*ldB+ j   ];
                tmp_j1 += A[i+l*ldA]*B[l*ldB+(j+1)];
                tmp_j2 += A[i+l*ldA]*B[l*ldB+(j+2)];
                tmp_j3 += A[i+l*ldA]*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i, ldA, B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }



}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_t_n(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*B[l+ j   *ldB];
                tmp_j1 += A[i*ldA+l]*B[l+(j+1)*ldB];
                tmp_j2 += A[i*ldA+l]*B[l+(j+2)*ldB];
                tmp_j3 += A[i*ldA+l]*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i*ldA, IndexType(1), B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }



}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_t_t(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*B[l*ldB+ j   ];
                tmp_j1 += A[i*ldA+l]*B[l*ldB+(j+1)];
                tmp_j2 += A[i*ldA+l]*B[l*ldB+(j+2)];
                tmp_j3 += A[i*ldA+l]*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i*ldA, IndexType(1), B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }



}


template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_c_n(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor_conj( k, A+i, ldA, &A_packed[ i*k ], 2*numElements);
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += conjugate(A[i+l*ldA])*B[l+ j   *ldB];
                tmp_j1 += conjugate(A[i+l*ldA])*B[l+(j+1)*ldB];
                tmp_j2 += conjugate(A[i+l*ldA])*B[l+(j+2)*ldB];
                tmp_j3 += conjugate(A[i+l*ldA])*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, A+i, ldA, B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_c_t(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor_conj( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += conjugate(A[i+l*ldA])*B[l*ldB+ j   ];
                tmp_j1 += conjugate(A[i+l*ldA])*B[l*ldB+(j+1)];
                tmp_j2 += conjugate(A[i+l*ldA])*B[l*ldB+(j+2)];
                tmp_j3 += conjugate(A[i+l*ldA])*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, A+i, ldA, B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }



}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_ct_n(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor_conj( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += conjugate(A[i*ldA+l])*B[l+ j   *ldB];
                tmp_j1 += conjugate(A[i*ldA+l])*B[l+(j+1)*ldB];
                tmp_j2 += conjugate(A[i*ldA+l])*B[l+(j+2)*ldB];
                tmp_j3 += conjugate(A[i*ldA+l])*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, A+i*ldA, IndexType(1), B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_ct_t(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor_conj( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += conjugate(A[i*ldA+l])*B[l*ldB+ j   ];
                tmp_j1 += conjugate(A[i*ldA+l])*B[l*ldB+(j+1)];
                tmp_j2 += conjugate(A[i*ldA+l])*B[l*ldB+(j+2)];
                tmp_j3 += conjugate(A[i*ldA+l])*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, A+i*ldA, IndexType(1), B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_n_c(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4_conj( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*conjugate(B[l+ j   *ldB]);
                tmp_j1 += A[i+l*ldA]*conjugate(B[l+(j+1)*ldB]);
                tmp_j2 += A[i+l*ldA]*conjugate(B[l+(j+2)*ldB]);
                tmp_j3 += A[i+l*ldA]*conjugate(B[l+(j+3)*ldB]);
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, B+k*ldB, IndexType(1), A+i, ldA, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_n_ct(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4_conj( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*conjugate(B[l*ldB+ j   ]);
                tmp_j1 += A[i+l*ldA]*conjugate(B[l*ldB+(j+1)]);
                tmp_j2 += A[i+l*ldA]*conjugate(B[l*ldB+(j+2)]);
                tmp_j3 += A[i+l*ldA]*conjugate(B[l*ldB+(j+3)]);
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, B+j, ldB, A+i, ldA, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_t_c(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4_conj( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*conjugate(B[l+ j   *ldB]);
                tmp_j1 += A[i*ldA+l]*conjugate(B[l+(j+1)*ldB]);
                tmp_j2 += A[i*ldA+l]*conjugate(B[l+(j+2)*ldB]);
                tmp_j3 += A[i*ldA+l]*conjugate(B[l+(j+3)*ldB]);
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, B+j*ldB, IndexType(1), A+i*ldA, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_t_ct(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4_conj( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*conjugate(B[l*ldB+ j   ]);
                tmp_j1 += A[i*ldA+l]*conjugate(B[l*ldB+(j+1)]);
                tmp_j2 += A[i*ldA+l]*conjugate(B[l*ldB+(j+2)]);
                tmp_j3 += A[i*ldA+l]*conjugate(B[l*ldB+(j+3)]);
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*tmp_j0;
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*tmp_j1;
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*tmp_j2;
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*tmp_j3;
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dot(k, B+j, ldB, A+i*ldA, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*tmp;
        }
    }



}


template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_c_c(IndexType m, IndexType n, IndexType k,
             const T &alpha, const MA *A, IndexType ldA,
             const MB *B, IndexType ldB,
             const T &beta,  T *C, IndexType ldC,
             MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4_conj( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor_conj( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*B[l+ j   *ldB];
                tmp_j1 += A[i+l*ldA]*B[l+(j+1)*ldB];
                tmp_j2 += A[i+l*ldA]*B[l+(j+2)*ldB];
                tmp_j3 += A[i+l*ldA]*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*conjugate(tmp_j0);
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*conjugate(tmp_j1);
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*conjugate(tmp_j2);
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*conjugate(tmp_j3);
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i, ldA, B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*conjugate(tmp);
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_c_ct(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4_conj( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {

                PackMatrixColToColMajor_conj( k, A+i, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i+l*ldA]*B[l*ldB+ j   ];
                tmp_j1 += A[i+l*ldA]*B[l*ldB+(j+1)];
                tmp_j2 += A[i+l*ldA]*B[l*ldB+(j+2)];
                tmp_j3 += A[i+l*ldA]*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*conjugate(tmp_j0);
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*conjugate(tmp_j1);
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*conjugate(tmp_j2);
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*conjugate(tmp_j3);
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i, ldA, B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*conjugate(tmp);
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_ct_c(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToRowMajor_4_conj( k, B+j*ldB, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor_conj( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*B[l+ j   *ldB];
                tmp_j1 += A[i*ldA+l]*B[l+(j+1)*ldB];
                tmp_j2 += A[i*ldA+l]*B[l+(j+2)*ldB];
                tmp_j3 += A[i*ldA+l]*B[l+(j+3)*ldB];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*conjugate(tmp_j0);
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*conjugate(tmp_j1);
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*conjugate(tmp_j2);
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*conjugate(tmp_j3);
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i*ldA, IndexType(1), B+j*ldB, IndexType(1), tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*conjugate(tmp);
        }
    }




}

template <typename IndexType, typename T, typename MA, typename MB>
void
sub_gemm_ct_ct(IndexType m, IndexType n, IndexType k,
              const T &alpha, const MA *A, IndexType ldA,
              const MB *B, IndexType ldB,
              const T &beta,  T *C, IndexType ldC,
              MA *A_packed, MB *B_packed, bool first_call)
{
    const int numElements = Intrinsics<T, DEFAULT_INTRINSIC_LEVEL>::numElements;
    IndexType i, j;



    for ( j=0; j+3<n; j+=4 ){

        PackMatrixColToColMajor_4_conj( k, B+j, ldB, &B_packed[ j*k ] );
        for ( i=0; i+2*numElements-1<m; i+=2*numElements ){

            if ( j == 0 ) {
                PackMatrixColToRowMajor_conj( k, A+i*ldA, ldA, &A_packed[ i*k ], 2*numElements );
            }
            kernel_gemm_2numElementsx4( k, alpha, &A_packed[ i*k ], 2*numElements, &B_packed[ j*k ], k, beta, C+i+j*ldC, ldC );
        }

        for (; i<m; ++i) {
            T tmp_j0 = 0,
              tmp_j1 = 0,
              tmp_j2 = 0,
              tmp_j3 = 0;
            for (IndexType l=0; l<k; ++l) {
                // C(i,j) = A(i,l)*B(l,j)
                tmp_j0 += A[i*ldA+l]*B[l*ldB+ j   ];
                tmp_j1 += A[i*ldA+l]*B[l*ldB+(j+1)];
                tmp_j2 += A[i*ldA+l]*B[l*ldB+(j+2)];
                tmp_j3 += A[i*ldA+l]*B[l*ldB+(j+3)];
            }
            C[i+ j   *ldC] = C[i+ j   *ldC] + alpha*conjugate(tmp_j0);
            C[i+(j+1)*ldC] = C[i+(j+1)*ldC] + alpha*conjugate(tmp_j1);
            C[i+(j+2)*ldC] = C[i+(j+2)*ldC] + alpha*conjugate(tmp_j2);
            C[i+(j+3)*ldC] = C[i+(j+3)*ldC] + alpha*conjugate(tmp_j3);
        }

    }

    for( ; j<n; ++j ){
        for (i=0; i<m; ++i) {
            T tmp(0);
            cxxblas::dotu(k, A+i*ldA, IndexType(1), B+j, ldB, tmp);
            C[i+j*ldC] = C[i+j*ldC] + alpha*conjugate(tmp);
        }
    }



}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_SUBGEMM_TCC
