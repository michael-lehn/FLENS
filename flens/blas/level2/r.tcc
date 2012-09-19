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

#ifndef FLENS_BLAS_LEVEL2_R_TCC
#define FLENS_BLAS_LEVEL2_R_TCC 1

#include <flens/blas/closures/debugclosure.h>
#include <flens/blas/level2/level2.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- GeneralMatrix, DenseVector ------------------------------------------------

//-- ger
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsGeMatrix<MA>::value,
         void>::Type
r(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    if ((A.numRows()==0) || (A.numCols()==0)) {
        A.resize(x.length(), y.length(), x.firstIndex(), y.firstIndex());
    }
    ASSERT(A.numRows()==x.length());
    ASSERT(A.numCols()==y.length());
    cxxblas::ger(A.order(),
                 A.numRows(), A.numCols(),
                 alpha,
                 x.data(), x.stride(),
                 y.data(), y.stride(),
                 A.data(), A.leadingDimension());
}

//-- geru
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsDenseVector<VY>::value
                 && IsGeMatrix<MA>::value,
         void>::Type
ru(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    r(alpha, x, y, A);
}

//-- gerc
template <typename ALPHA, typename VX, typename VY, typename MA>
    typename RestrictTo<IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value
                     && IsGeMatrix<MA>::value,
             void>::Type
rc(const ALPHA &alpha, const VX &x, const VY &y, MA &&A)
{
    if ((A.numRows()==0) || (A.numCols()==0)) {
        A.resize(x.length(), y.length(), x.firstIndex(), y.firstIndex());
    }
    ASSERT(A.numRows()==x.length());
    ASSERT(A.numCols()==y.length());
    cxxblas::gerc(A.order(),
                  A.numRows(), A.numCols(),
                  alpha,
                  x.data(), x.stride(),
                  y.data(), y.stride(),
                  A.data(), A.leadingDimension());
}

//-- HermitianMatrix, DenseVector ----------------------------------------------

//-- her
template <typename ALPHA, typename VX, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsHeMatrix<MA>::value,
         void>::Type
r(const ALPHA &alpha, const VX &x, MA &&A)
{
    if (A.dim()==0) {
        A.resize(x.length(), x.firstIndex());
    }
    ASSERT(A.dim()==x.length());
    cxxblas::her(A.order(),
                 A.upLo(),
                 A.dim(),
                 alpha,
                 x.data(), x.stride(),
                 A.data(), A.leadingDimension());
}

//-- hpr
template <typename ALPHA, typename VX, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsHpMatrix<MA>::value,
         void>::Type
r(const ALPHA &alpha, const VX &x, MA &&A)
{
    if (A.dim()==0) {
        A.resize(x.length(), x.firstIndex());
    }
    ASSERT(A.dim()==x.length());
    cxxblas::hpr(A.order(),
                 A.upLo(),
                 A.dim(),
                 alpha,
                 x.data(), x.stride(),
                 A.data());
}

//-- SymmetricMatrix, DenseVector ----------------------------------------------

//-- spr
template <typename ALPHA, typename VX, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsSpMatrix<MA>::value,
         void>::Type
r(const ALPHA &alpha, const VX &x, MA &&A)
{
    if (A.dim()==0) {
        A.resize(x.length(), x.firstIndex());
    }
    ASSERT(A.dim()==x.length());
    cxxblas::spr(A.order(),
                 A.upLo(),
                 A.dim(),
                 alpha,
                 x.data(), x.stride(),
                 A.data());
}

//-- syr
template <typename ALPHA, typename VX, typename MA>
typename RestrictTo<IsDenseVector<VX>::value
                 && IsSyMatrix<MA>::value,
         void>::Type
r(const ALPHA &alpha, const VX &x, MA &&A)
{
    if (A.dim()==0) {
        A.resize(x.length(), x.firstIndex());
    }
    ASSERT(A.dim()==x.length());
    cxxblas::syr(A.order(),
                 A.upLo(),
                 A.dim(),
                 alpha,
                 x.data(), x.stride(),
                 A.data(), A.leadingDimension());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_R_TCC
