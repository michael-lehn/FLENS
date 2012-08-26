/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef CXXBLAS_SPARSELEVEL2_HECCSMV_TCC
#define CXXBLAS_SPARSELEVEL2_HECCSMV_TCC 1

#include <cxxblas/auxiliary/auxiliary.h>
#include <cxxblas/typedefs.h>

#define HAVE_CXXBLAS_HECRSMV 1

namespace cxxblas {

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
heccsmv(StorageUpLo      upLo,
        IndexType        n,
        const ALPHA      &alpha,
        const MA         *A,
        const IndexType  *ia,
        const IndexType  *ja,
        const VX         *x,
        const BETA       &beta,
        VY               *y)
{
    using cxxblas::conjugate;

//
//  The correct index base of the CCS matrix is stored in first Element of ja
//
    --ja;
    ia -= ja[1];
    A  -= ja[1];

//
//  Let x, y be one-based; x_, y_ get correct index base ja[1]
//
    const VX *x_ = x - ja[1];
    VY       *y_ = y - ja[1];
    --x;
    --y;

    if (beta==BETA(0)) {
        for (int i=1; i<=n; ++i) {
            y[i] = 0;
        }
    } else if (beta!=BETA(1)) {
        for (int i=1; i<=n; ++i) {
            y[i] *= beta;
        }
    }

    if (upLo==Lower) {
        for (int j=1, J=ja[1]; j<=n; ++j, ++J) {
            if (ja[j]<ja[j+1]) {
                int k=ja[j];
                if (ia[k]==J) {
                    y[j]      += alpha*cxxblas::real(A[k])*x_[ia[k]];
                } else {
                    y[j]      += alpha*conjugate(A[k])*x_[ia[k]];
                    y_[ia[k]] += alpha*A[k]*x[j];
                }
                for (k=ja[j]+1; k<ja[j+1]; ++k) {
                    y[j]      += alpha*conjugate(A[k])*x_[ia[k]];
                    y_[ia[k]] += alpha*A[k]*x[j];
                }
            }
        }
    } else {
        for (int j=1, J=ja[1]; j<=n; ++j, ++J) {
            if (ja[j]<ja[j+1]) {
                int k;
                for (k=ja[j]; k<ja[j+1]-1; ++k) {
                    y[j]      += alpha*conjugate(A[k])*x_[ia[k]];
                    y_[ia[k]] += alpha*A[k]*x[j];
                }
                if (ia[k]==J) {
                    y[j]      += alpha*cxxblas::real(A[k])*x_[ia[k]];
                } else {
                    y[j]      += alpha*conjugate(A[k])*x_[ia[k]];
                    y_[ia[k]] += alpha*A[k]*x[j];
                }
            }
        }
    }
}

} // namespace cxxblas

#endif // CXXBLAS_SPARSELEVEL2_HECCSMV_TCC
