/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef CXXBLAS_LEVEL2_TPSV_TCC
#define CXXBLAS_LEVEL2_TPSV_TCC 1

namespace cxxblas {

template <typename IndexType, typename MA, typename VX>
void
tpsv_generic(StorageOrder order, StorageUpLo upLo,
             Transpose transA, Diag diag,
             IndexType n,
             const MA *A,
             VX *x, IndexType incX)
{
    if (order==ColMajor) {
        transA = Transpose(transA^Trans);
        upLo = (upLo==Upper) ? Lower : Upper;
        tpsv_generic(RowMajor, upLo, transA, diag, n, A, x, incX);
        return;
    }
    if (transA==NoTrans) {
        if (upLo==Upper) {
            if (diag==NonUnit) {
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    VX _x;
                    dotu_generic(n-1-i, A+i*(2*n-i+1)/2+1, IndexType(1),
                                        x+iX+incX, incX, _x);
                    x[iX] -= _x;
                    x[iX] /= *(A+i*(2*n-i+1)/2);
                }
            } else { /* diag==Unit */
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    VX _x;
                    dotu_generic(n-1-i, A+i*(2*n-i+1)/2+1, IndexType(1),
                                        x+iX+incX, incX, _x);
                    x[iX] -= _x;
                }
            }
        } else { /* upLo==Lower */
            if (diag==NonUnit) {
                for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
                    VX _x;
                    dotu_generic(i, A+i*(i+1)/2, IndexType(1),
                                    x, incX, _x);
                    x[iX] -= _x;
                    x[iX] /= *(A+i*(i+3)/2);
                }
            } else { /* diag==Unit */
                for (IndexType i=1, iX=i*incX; i<n; ++i, iX+=incX) {
                    VX _x;
                    dotu_generic(i, A+i*(i+1)/2, IndexType(1),
                                    x, incX, _x);
                    x[iX] -= _x;
                }
            }
        }
    }
    if (transA==Conj) {
        if (upLo==Upper) {
            if (diag==NonUnit) {
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    VX _x;
                    dot_generic(n-1-i, A+i*(2*n-i+1)/2+1, IndexType(1),
                                       x+iX+incX, incX, _x);
                    x[iX] -= _x;
                    x[iX] /= conjugate(A[i*(2*n-i+1)/2]);
                }
            } else { /* diag==Unit */
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    VX _x;
                    dot_generic(n-1-i, A+i*(2*n-i+1)/2+1, IndexType(1),
                                       x+iX+incX, incX, _x);
                    x[iX] -= _x;
                }
            }
        } else { /* upLo==Lower */
            if (diag==NonUnit) {
                for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
                    VX _x;
                    dot_generic(i, A+i*(i+1)/2, IndexType(1),
                                   x, incX, _x);
                    x[iX] -= _x;
                    x[iX] /= conjugate(A[i*(i+3)/2]);
                }
            } else { /* diag==Unit */
                for (IndexType i=1, iX=i*incX; i<n; ++i, iX+=incX) {
                    VX _x;
                    dot_generic(i, A+i*(i+1)/2, IndexType(1),
                                   x, incX, _x);
                    x[iX] -= _x;
                }
            }
        }
    }
    if (transA==Trans) {
        if (upLo==Upper) {
            if (diag==NonUnit) {
                for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
                    x[iX] /= *(A+i*(2*n-i+1)/2);
                    axpy_generic(n-1-i, -x[iX], A+i*(2*n-i+1)/2+1, IndexType(1),
                                                x+iX+incX, incX);
                }
            } else { /* diag==Unit */
                for (IndexType i=0, iX=0; i<n-1; ++i, iX+=incX) {
                    axpy_generic(n-1-i, -x[iX], A+i*(2*n-i+1)/2+1, IndexType(1),
                                                x+iX+incX, incX);
                }
            }
        } else { /* upLo==Lower */
            if (diag==NonUnit) {
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    x[iX] /= *(A+i*(i+3)/2);
                    axpy_generic(i, -x[iX], A+i*(i+1)/2, IndexType(1),
                                            x, incX);
                }
            } else {
                for (IndexType i=n-1, iX=i*incX; i>0; --i, iX-=incX) {
                    axpy_generic(i, -x[iX], A+i*(i+1)/2, IndexType(1),
                                            x, incX);
                }
            }
        }
    }
    if (transA==ConjTrans) {
        if (upLo==Upper) {
            if (diag==NonUnit) {
                for (IndexType i=0, iX=0; i<n; ++i, iX+=incX) {
                    x[iX] /= conjugate(A[i*(2*n-i+1)/2]);
                    acxpy_generic(n-1-i, -x[iX], A+i*(2*n-i+1)/2+1, IndexType(1),
                                                 x+iX+incX, incX);
                }
            } else { /* diag==Unit */
                for (IndexType i=0, iX=0; i<n-1; ++i, iX+=incX) {
                    acxpy_generic(n-1-i, -x[iX], A+i*(2*n-i+1)/2+1, IndexType(1),
                                                 x+iX+incX, incX);
                }
            }
        } else { /* upLo==Lower */
            if (diag==NonUnit) {
                for (IndexType i=n-1, iX=i*incX; i>=0; --i, iX-=incX) {
                    x[iX] /= conjugate(A[i*(i+3)/2]);
                    acxpy_generic(i, -x[iX], A+i*(i+1)/2, IndexType(1),
                                             x, incX);
                }
            } else {
                for (IndexType i=n-1, iX=i*incX; i>0; --i, iX-=incX) {
                    acxpy_generic(i, -x[iX], A+i*(i+1)/2, IndexType(1),
                                             x, incX);
                }
            }
        }
    }
}

template <typename IndexType, typename MA, typename VX>
void
tpsv(StorageOrder order, StorageUpLo upLo,
     Transpose transA, Diag diag,
     IndexType n,
     const MA *A,
     VX *x, IndexType incX)
{
    CXXBLAS_DEBUG_OUT("tpsv_generic");

    if (incX<0) {
        x -= incX*(n-1);
    }
    tpsv_generic(order, upLo, transA, diag, n, A, x, incX);
}

} // namespace flens

#endif // CXXBLAS_LEVEL2_TPSV_TCC

