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
 *
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
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

//== generic lapack implementation =============================================

template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg_generic(IndexType n, ALPHA &alpha, DenseVector<VX> &x, TAU &tau)
{
    using std::abs;

    typedef typename DenseVector<VX>::ElementType   T;

    if (n<=1) {
        tau = TAU(0);
        return;
    }

    T xNorm = blas::nrm2(x);
    if (xNorm==T(0)) {
//
//      H  =  I
//
        tau = TAU(0);
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

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK

template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg_native(IndexType n, ALPHA &alpha, DenseVector<VX> &x, TAU &tau)
{
    typedef typename DenseVector<VX>::ElementType  T;

    INTEGER N       = n;
    INTEGER INCX    = x.inc();

    if (IsSame<T, DOUBLE>::value) {
        LAPACK_IMPL(dlarfg)(&N,
                            &alpha,
                            x.data(),
                            &INCX,
                            &tau);
    } else {
        ASSERT(0);
    }
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg(IndexType n, ALPHA &alpha, DenseVector<VX> &x, TAU &tau)
{
    LAPACK_DEBUG_OUT("larfg");

//
//  Test the input parameters
//
    ASSERT(x.firstIndex()==1);
    ASSERT(x.inc()>0);
    ASSERT(x.length()<=n);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    ALPHA                               _alpha  = alpha;
    typename DenseVector<VX>::NoView    _x      = x;
    TAU                                 _tau    = tau;
#   endif

//
//  Call implementation
//
    larfg_generic(n, alpha, x, tau);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    larfg_native(n, _alpha, _x, _tau);

    bool failed = false;
    if (alpha!=_alpha) {
        std::cerr << "CXXLAPACK:  alpha = " << alpha << std::endl;
        std::cerr << "F77LAPACK: _alpha = " << _alpha << std::endl;
        std::cerr << "CXXLAPACK:  alpha - _alpha= "
                  << alpha - _alpha << std::endl;
        failed = true;
    }

    if (! isIdentical(x, _x, " x", "x_")) {
        std::cerr << "CXXLAPACK:  x = " << x << std::endl;
        std::cerr << "F77LAPACK: _x = " << _x << std::endl;
        failed = true;
    }

    if (tau!=_tau) {
        std::cerr << "CXXLAPACK:  tau = " << tau << std::endl;
        std::cerr << "F77LAPACK: _tau = " << _tau << std::endl;
        std::cerr << "CXXLAPACK:  tau - _tau= "
                  << tau - _tau << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename ALPHA, typename VX, typename TAU>
void
larfg(IndexType n, ALPHA &alpha, VX &&x, TAU &tau)
{
    larfg(n, alpha, x, tau);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LARFG_TCC
