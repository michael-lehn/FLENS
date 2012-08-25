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

#ifndef CXXBLAS_SPARSELEVEL2_GECRSMV_TCC
#define CXXBLAS_SPARSELEVEL2_GECRSMV_TCC 1

#include <cxxblas/typedefs.h>

#define HAVE_CXXBLAS_GECRSMV 1

namespace cxxblas {

template <typename IndexType, typename ALPHA, typename MA, typename VX,
          typename BETA, typename VY>
void
gecrsmv(Transpose        trans,
        IndexType        m,
        IndexType        n,
        const ALPHA      &alpha,
        const MA         *A,
        const IndexType  *ia,
        const IndexType  *ja,
        const VX         *x,
        const BETA       &beta,
        VY               *y)
{
    const bool init  = (beta==BETA(0));
    const bool scale = (beta!=BETA(0) && beta!=BETA(1));

    ASSERT(incX==IndexType(1));
    ASSERT(incY==IndexType(1));
    /*
    // If matrix is Zero based:
    // set index base of y to 1
    --y;
    --ia;
    */

    /*
    // If matrix is One based
    --y;
    --ia;
    --ja;
    --x;
    --A;
    */

//
//  Index base of the CRS matrix is stored in first Element of ia
//
    --ia;
    ja -= ia[1];
    A  -= ia[1];

    if (trans==NoTrans) {
//
//      Make y one-based; set correct index base for x
//
        --y;
        x  -= ia[1];

        if (init) {
            for (int i=1; i<=m; ++i) {
                y[i] = VY(0);
                for (int k=ia[i]; k<ia[i+1]; ++k) {
                    y[i] += alpha*A[k]*x[ja[k]];
                }
            }
        } else if (scale) {
            for (int i=1; i<=m; ++i) {
                y[i] *= beta;
                for (int k=ia[i]; k<ia[i+1]; ++k) {
                    y[i] += alpha*A[k]*x[ja[k]];
                }
            }
        } else {
            for (int i=1; i<=m; ++i) {
                for (int k=ia[i]; k<ia[i+1]; ++k) {
                    y[i] += alpha*A[k]*x[ja[k]];
                }
            }
        }
    } else {
//
//      Make y one-based
//
        --y;
        if (init) {
            for (int i=1; i<=n; ++i) {
                y[i] = VY(0);
            }
        } else if (scale) {
            for (int i=1; i<=n; ++i) {
                y[i] *=beta;
            }
        }
//
//      Set corret index base for y; make x one-based
//
        y += 1-ia[1];
        --x;

        for (int i=1; i<=m; ++i) {
            for (int k=ia[i]; k<ia[i+1]; ++k) {
                y[ja[k]] += alpha*A[k]*x[i];
            }
        }
    }
}

} // namespace cxxblas

#endif // CXXBLAS_SPARSELEVEL2_GECRSMV_TCC
