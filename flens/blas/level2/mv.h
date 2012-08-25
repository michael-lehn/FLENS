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

#ifndef FLENS_BLAS_LEVEL2_MV_H
#define FLENS_BLAS_LEVEL2_MV_H 1

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//== GeneralMatrix - Vector products ===========================================

//-- gemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsGeMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
       const BETA &beta, VY &&y);

//-- gecrsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsGeCRSMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
       const BETA &beta, VY &&y);

//-- geccsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsGeCCSMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(Transpose trans, const ALPHA &alpha, const MA &A, const VX &x,
       const BETA &beta, VY &&y);

//== TriangularMatrix - Vector products ========================================

//-- trmv
template <typename MA, typename VX>
    typename RestrictTo<IsTrMatrix<MA>::value
                     && IsDenseVector<VX>::value,
             void>::Type
    mv(Transpose trans, const MA &A, VX &&x);


//== SymmetricMatrix - Vector products =========================================

//-- symv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsSyMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

//-- sycrsmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsSyCRSMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);


//== HermitianMatrix - Vector products =========================================

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsHeMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VY>::value,
             void>::Type
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_H
