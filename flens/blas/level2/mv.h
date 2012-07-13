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
    void
    mv(Transpose trans,
       const ALPHA &alpha, const GeMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(Transpose trans,
       const ALPHA &alpha, const GbMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

//== TriangularMatrix - Vector products ========================================

//-- tbmv
template <typename MA, typename VX>
    void
    mv(Transpose trans, const TbMatrix<MA> &A, DenseVector<VX> &x);

//-- trmv
template <typename MA, typename VX>
    void
    mv(Transpose trans, const TrMatrix<MA> &A, DenseVector<VX> &x);

//-- tpmv
template <typename MA, typename VX>
    void
    mv(Transpose trans, const TpMatrix<MA> &A, DenseVector<VX> &x);
    
//== SymmetricMatrix - Vector products =========================================

//-- sbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const SbMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

//-- spmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const SpMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);
    
//-- symv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const SyMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

//== HermitianMatrix - Vector products =========================================

//-- hbmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const HbMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const HeMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);
    
//-- hpmv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const HpMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);

//== forwarding ================================================================

//-- GeneralMatrix - Vector products
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsGeneralMatrix<MA>::value &&
                       !IsClosure<MA>::value && !IsClosure<VX>::value &&
                        IsSame<VY, typename VY::Impl>::value,
             void>::Type
    mv(Transpose trans,
       const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

//-- TriangularMatrix - Vector products
template <typename MA, typename VX>
    typename RestrictTo<IsTriangularMatrix<MA>::value &&
                       !IsClosure<MA>::value &&
                        IsSame<VX, typename VX::Impl>::value,
             void>::Type
    mv(Transpose trans, const MA &A, VX &&x);

//-- Symmetric Matrix - Vector products
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsSymmetricMatrix<MA>::value &&
                       !IsClosure<MA>::value && !IsClosure<VX>::value &&
                        IsSame<VY, typename VY::Impl>::value,
             void>::Type
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

//-- Hermitian Matrix - Vector products
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    typename RestrictTo<IsHermitianMatrix<MA>::value &&
                       !IsClosure<MA>::value && !IsClosure<VX>::value &&
                        IsSame<VY, typename VY::Impl>::value,
             void>::Type
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);



} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_H
