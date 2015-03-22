/*
 *   Copyright (c) 2012, 2015 Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL2_R2_TCC
#define FLENS_BLAS_LEVEL2_R2_TCC 1

#include <flens/blas/closures/closures.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- HermitianMatrix, DenseVector ----------------------------------------------

//-- her2
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsHeMatrix<MA>::value,
         void>::Type
r2(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    const bool lowerA    = (A.upLo()==Lower);

    ASSERT(x.length()==y.length());

    if (A.dim()==0) {
        ASSERT(x.firstIndex()==y.firstIndex());
        A.resize(x.length(), x.firstIndex());
    }

    ASSERT(A.dim()==x.length());

    cxxblas::her2(A.dim(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  lowerA, A.data(), A.strideRow(), A.strideCol());
}

//-- hpr2
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsHpMatrix<MA>::value,
         void>::Type
r2(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    ASSERT(x.length()==y.length());

    if (A.dim()==0) {
        ASSERT(x.firstIndex()==y.firstIndex());
        A.resize(x.length(), x.firstIndex());
    }

    ASSERT(A.dim()==x.length());

    cxxblas::hpr2(A.dim(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  colMajorA, lowerA, A.data());
}

//-- SymmetricMatrix, DenseVector ----------------------------------------------

//-- spr2
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsSpMatrix<MA>::value,
         void>::Type
r2(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    const bool colMajorA = (A.order()==ColMajor);
    const bool lowerA    = (A.upLo()==Lower);

    ASSERT(x.length()==y.length());

    if (A.dim()==0) {
        ASSERT(x.firstIndex()==y.firstIndex());
        A.resize(x.length(), x.firstIndex());
    }

    ASSERT(A.dim()==x.length());

    cxxblas::spr2(A.dim(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  colMajorA, lowerA, A.data());
}

//-- syr2
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsSyMatrix<MA>::value,
         void>::Type
r2(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    const bool lowerA    = (A.upLo()==Lower);

    ASSERT(x.length()==y.length());

    if (A.dim()==0) {
        ASSERT(x.firstIndex()==y.firstIndex());
        A.resize(x.length(), x.firstIndex());
    }

    ASSERT(A.dim()==x.length());

    cxxblas::syr2(A.dim(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  lowerA, A.data(), A.strideRow(), A.strideCol());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_R2_TCC
