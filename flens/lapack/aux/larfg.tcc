/*
 *   Copyright (c) 2011, Michael Lehn
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

/* Baesed on
      SUBROUTINE xLARFG( N, ALPHA, X, INCX, TAU )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_AUX_LARFG_TCC
#define FLENS_LAPACK_AUX_LARFG_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

using std::abs;

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg(IndexType n, ALPHA &alpha, VX &&x, TAU &tau)
{
    larfg(n, alpha, x, tau);
}

//-- larfg ---------------------------------------------------------------------
template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg(IndexType n, ALPHA &alpha, DenseVector<VX> &x, TAU &tau)
{
    ASSERT(x.inc()>0);
    ASSERT(x.length()<=n);
    ASSERT(x.firstIndex()==1);

    typedef typename DenseVector<VX>::ElementType   T;

    if (n<=1) {
        tau = 0;
        return;
    }

    T xNorm = blas::nrm2(x);
    if (xNorm==T(0)) {
//
//      H  =  I
//
        tau = 0;
    } else {
//
//      general case
//
        T beta = -sign(lapy2(alpha, xNorm), alpha);
        T safeMin = lamch<T>(SafeMin) / lamch<T>(Eps);

        IndexType count=0;
        if (abs(beta)<safeMin) {
//
//          XNORM, BETA may be inaccurate; scale X and recompute them
//
            T rSafeMin = T(1)/safeMin;
            do {
                ++count;
                blas::scal(rSafeMin, x);
                beta *= rSafeMin;
                alpha *= rSafeMin;
            } while (abs(beta)<safeMin);
//
//          New BETA is at most 1, at least SAFMIN
//
            xNorm = blas::nrm2(x);
            beta = -sign(lapy2(alpha, xNorm), alpha);
        }
        tau = (beta-alpha) / beta;
        blas::scal(T(1)/(alpha-beta), x);
//
//      If ALPHA is subnormal, it may lose relative accuracy
//
        for (IndexType j=1; j<=count; ++j) {
            beta *= safeMin;
        }
        alpha = beta;
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFG_TCC
