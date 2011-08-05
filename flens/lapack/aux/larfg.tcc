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

/*
 *   Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 *   
 *   $COPYRIGHT$
 *   
 *   Additional copyrights may follow
 *   
 *   $HEADER$
 *   
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are
 *   met:
 *   
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer. 
 *     
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer listed
 *     in this license in the documentation and/or other materials
 *     provided with the distribution.
 *     
 *   - Neither the name of the copyright holders nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
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

        // hack: gmp returns safeMin==0
        T eps = lamch<T>(Eps);
        T safeMin = lamch<T>(SafeMin);
        if (safeMin>0) {
            safeMin /= eps;
        }

        if (abs(beta)<safeMin) {
//
//          xNorm, beta may be inaccurate; scale x and recompute them
//
            const T rSafeMin = T(1)/safeMin;

            IndexType count = 0;
            do {
                ++count;
                blas::scal(rSafeMin, x);
                beta = beta*rSafeMin;
                alpha = alpha*rSafeMin;
                
            } while (abs(beta) < safeMin);
//
//          New beta is at most 1, at least safeMin
//
            xNorm = blas::nrm2(x);
            beta = -sign(lapy2(alpha, xNorm), alpha);
            tau = (beta-alpha) / beta;
            blas::scal(1/(alpha-beta), x);
//
//          If alpha is subnormal, it may lose relative accuracy
//
            alpha = beta;
            for (IndexType j=1; j<=count; ++j) {
                alpha *= safeMin;
            }
        } else {
            tau = (beta-alpha) / beta;
            blas::scal(1/(alpha-beta), x);
            alpha = beta;
        }
    }
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFG_TCC
