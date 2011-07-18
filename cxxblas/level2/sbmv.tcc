/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef CXXBLAS_LEVEL2_SBMV_TCC
#define CXXBLAS_LEVEL2_SBMV_TCC 1

#include <complex>
#include <cxxblas/level1/level1.h>

using std::min;
using std::max;

namespace cxxblas {

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
sbmv_generic(StorageOrder order, StorageUpLo upLo,
             IndexType n, IndexType k,
             const ALPHA &alpha,
             const MA *A, IndexType ldA,
             const VX *x, IndexType incX,
             const BETA &beta,
             VY *y, IndexType incY)
{
    if (order==ColMajor) {
        upLo = (upLo==Upper) ? Lower : Upper;
        sbmv_generic(RowMajor, upLo, n, k, alpha, A, ldA,
                     x, incX, beta, y, incY);
        return;
    }
    scal_generic(n, beta, y, incY);
    if (upLo==Upper) {
        for (IndexType i=0, iX=0, iY=0; i<n; ++i, iX+=incX, iY+=incY) {
            IndexType len = min(k+1, n-i);

            VY _y;
            dot_generic(len, A+ldA*i, IndexType(1),
                             x+iX, IndexType(incX),
                             _y);
            y[iY] += alpha*_y;

            axpy_generic(len-1, x[iX] * alpha,
                                A+ldA*i+1, IndexType(1),
                                y+iY+incY, incY);
        }
    } else {
        for (IndexType i=0, iY=0; i<n; ++i, iY+=incY) {
            IndexType iA = max(k-i, IndexType(0));
            IndexType len = min(k, i) + 1;
            IndexType _i = max(i-k, IndexType(0));

            VY _y;
            dot_generic(len, A+ldA*i+iA, IndexType(1),
                             x+_i*incX, IndexType(incX),
                             _y);
            y[iY] += alpha*_y;

            axpy_generic(len-1, x[i*incX] * alpha,
                                A+ldA*i+iA, IndexType(1),
                                y+_i*incY, incY);
        }
    }
}

//------------------------------------------------------------------------------

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
sbmv(StorageOrder order, StorageUpLo upLo,
     IndexType n, IndexType k,
     const ALPHA &alpha,
     const MA *A, IndexType ldA,
     const VX *x, IndexType incX,
     const BETA &beta,
     VY *y, IndexType incY)
{
    if (n==0) {
        return;
    }
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    sbmv_generic(order, upLo, n, k, alpha, A, ldA, x, incX, beta, y, incY);
}

} // namespace cxxblas

#endif // CXXBLAS_LEVEL2_SBMV_TCC
