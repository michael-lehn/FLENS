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
 *
 *   Algorithm based on the tutorial for fast matrix-matrix mulitplication:
 *
 *   http://wiki.cs.utexas.edu/rvdg/HowToOptimizeGemm
 *
 *   It is based on the papers:
 *
 *    Anatomy of high-performance matrix multiplication.
 *    Kazushige Goto, Robert A. van de Geijn.
 *    ACM Transactions on Mathematical Software (TOMS), 2008.
 *
 *
 *    BLIS: A Modern Alternative to the BLAS.
 *    Field G. Van Zee, Robert A. van de Geijn.
 *    ACM Transactions on Mathematical Software, submitted.
 *
 */

#ifndef PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC
#define PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC 1

#include <playground/cxxblas/intrinsics/level3/gemm/packmatrix.tcc>
#include <playground/cxxblas/intrinsics/level3/gemm/subgemm.tcc>
#include <playground/cxxblas/intrinsics/level3/gemm/kernelgemm.tcc>

namespace cxxblas {

#ifdef USE_INTRINSIC

template <typename IndexType, typename T, typename MA, typename MB>
typename flens::RestrictTo< flens::IsIntrinsicsCompatible<T>::value &&
                            flens::IsIntrinsicsCompatible<MA>::value &&
                            flens::IsIntrinsicsCompatible<MB>::value &&
                            (flens::IsComplex<T>::value ||
                                (flens::IsReal<MA>::value &&
                                 flens::IsReal<MB>::value) ),
                            void>::Type
gemm(StorageOrder order, Transpose transA, Transpose transB,
     IndexType m, IndexType n, IndexType k,
     const T &alpha,
     const MA *A, IndexType ldA,
     const MB *B, IndexType ldB,
     const T &beta,
     T *C, IndexType ldC)
{
    CXXBLAS_DEBUG_OUT("gemm_intrinsics [" INTRINSIC_NAME "]");

    using std::min;

    if ((m==0) || (n==0)) {
        return;
    }

    if (order==RowMajor) {
        gemm(ColMajor, transB, transA,
                     n, m, k, alpha,
                     B, ldB, A, ldA,
                     beta,
                     C, ldC);
        return;
    }

    cxxblas::gescal(cxxblas::StorageOrder::ColMajor, m, n, beta, C, ldC);

    IndexType kc = BLOCKSIZE_GEMM_K;
    IndexType mc = BLOCKSIZE_GEMM_M;

    // Allocate memory
    MA *packedA = new MA[ min(k,kc) * min(m,mc) ];
    MB *packedB = new MB[ min(k,kc) * n ];

    if ( transA==NoTrans ) {

        if ( transB==NoTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_n_n(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==Conj ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_n_c(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0 );
                }
            }

        } else if ( transB==Trans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_n_t(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p*ldB, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==ConjTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_n_ct(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                  B+p*ldB, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }

        }

    } else if ( transA==Conj ) {

        if ( transB==NoTrans || transB==Conj ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_c_n(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==Conj ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_c_c(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==Trans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_c_t(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                 B+p*ldB, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==ConjTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_c_ct(ib, n, pb, alpha, A+i+p*ldA, ldA,
                                  B+p*ldB, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }

        }

    } else if ( transA==Trans ) {

        if ( transB==NoTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_t_n(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==Conj ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_t_c(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                 B+p, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }

        } else if ( transB==Trans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_t_t(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                 B+p*ldB, ldB, T(1), C+i, ldC,
                                 packedA, packedB, i==0);
                }
            }
        } else if ( transB==ConjTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_t_ct(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                  B+p*ldB, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }
        }

    } else if ( transA==ConjTrans ) {

        if ( transB==NoTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_ct_n(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                  B+p, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }

        } else if ( transB==Conj ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_ct_c(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                  B+p, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }

        } else if ( transB==Trans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_ct_t(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                  B+p*ldB, ldB, T(1), C+i, ldC,
                                  packedA, packedB, i==0);
                }
            }
        } else if ( transB==ConjTrans ) {

            for ( IndexType p=0; p<k; p+=kc ){
                IndexType pb = min( k-p, kc );

                for (IndexType i=0; i<m; i+=mc ){
                    IndexType ib = min( m-i, mc );
                    sub_gemm_ct_ct(ib, n, pb, alpha, A+p+i*ldA, ldA,
                                   B+p*ldB, ldB, T(1), C+i, ldC,
                                   packedA, packedB, i==0);
                }
            }
        }
    }

    // free memory again
    delete[] packedA;
    delete[] packedB;

}

#endif // USE_INTRINSIC

} // namespace cxxblas

#endif // PLAYGROUND_CXXBLAS_INTRINSICS_LEVEL3_GEMM_TCC
