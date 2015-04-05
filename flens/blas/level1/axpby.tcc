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

#ifndef FLENS_BLAS_LEVEL1_AXPBY_TCC
#define FLENS_BLAS_LEVEL1_AXPBY_TCC 1

#include <flens/matrixtypes/matrixtypes.h>
#include <flens/storage/storage.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- BLAS Level 1 extensions ---------------------------------------------------

//-- axpby
template <typename ALPHA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
axpby(const ALPHA &alpha, const VX &x, const BETA &beta, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_AXPBY(alpha, x, beta, y);

    if (y.length()==0) {
//
//      So we allow  y = beta*y + alpha*x  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

    cxxblas::axpby(x.length(), alpha,
                   x.data(), x.stride(),
                   beta,
                   y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

//-- acxpby
template <typename ALPHA, typename VX, typename BETA, typename VY>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value,
         void>::Type
acxpby(const ALPHA &alpha, const VX &x, const BETA &beta, VY &&y)
{
    FLENS_BLASLOG_SETTAG("--> ");
    FLENS_BLASLOG_BEGIN_AXPY(alpha, x, y);

    if (y.length()==0) {
//
//      So we allow  y = beta*y + alpha*x  for an empty vector y
//
        typedef typename RemoveRef<VY>::Type   VectorY;
        typedef typename VectorY::ElementType  T;
        const T  Zero(0);

        y.resize(x, Zero);
    }
    ASSERT(y.length()==x.length());

    cxxblas::acxpby(x.length(), alpha,
                    x.data(), x.stride(),
                    beta,
                    y.data(), y.stride());

    FLENS_BLASLOG_END;
    FLENS_BLASLOG_UNSETTAG;
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL1_AXPBY_TCC
