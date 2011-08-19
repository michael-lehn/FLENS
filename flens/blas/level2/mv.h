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
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- forwarding: GeneralMatrix - Vector products -------------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

//-- product type: GeneralMatrix - Vector products -----------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const GeneralMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

//-- gemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const GeMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);


//-- forwarding: Hermitian Matrix - Vector products ----------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const MA &A, const VX &x, const BETA &beta, VY &&y);

//-- product type: HermitianMatrix - Vector products ---------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const HermitianMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

//-- hemv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const HeMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);


//-- forwarding: SymmetricMatrix - Vector products -----------------------------
// -> is identical with forwarding of Hermitian Matrix - Vector products

//-- product type: SymmetricMatrix - Vector products ---------------------------
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const SymmetricMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

//-- symv
template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(const ALPHA &alpha, const SyMatrix<MA> &A, const DenseVector<VX> &x,
       const BETA &beta, DenseVector<VY> &y);


//-- forwarding: TriangularMatrix - Vector products ----------------------------
template <typename MA, typename VX>
    void
    mv(cxxblas::Transpose trans,  const MA &A, VX &&x);

//-- product type: TriangularMatrix - Vector products --------------------------
template <typename MA, typename VX>
    void
    mv(cxxblas::Transpose trans,  const TriangularMatrix<MA> &A, Vector<VX> &x);

//-- trmv
template <typename MA, typename VX>
    void
    mv(cxxblas::Transpose trans, const TrMatrix<MA> &A, DenseVector<VX> &x);


} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL2_MV_H
