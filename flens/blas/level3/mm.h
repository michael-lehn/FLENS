/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_MM_H
#define FLENS_BLAS_LEVEL3_MM_H

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//== product type: GeneralMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB,
       const ALPHA &alpha,
       const MA &A, const MB &B,
       const BETA &beta,
       MC &&C);

//-- common interface ----------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB,
       const ALPHA &alpha,
       const GeneralMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta,
       GeneralMatrix<MC> &C);

//-- gemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB,
       const ALPHA &alpha,
       const GeMatrix<MA> &A, const GeMatrix<MB> &B,
       const BETA &beta,
       GeMatrix<MC> &C);


//== product type: HermitianMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Side side, const ALPHA &alpha, const MA &A, const MB &B,
       const BETA &beta, MC &&C);

//-- common interface ----------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Side side,
       const ALPHA &alpha,
       const HermitianMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, GeneralMatrix<MC> &C);

//-- hemm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Side side,
       const ALPHA &alpha, const HeMatrix<MA> &A, const GeMatrix<MB> &B,
       const BETA &beta, GeMatrix<MC> &C);


//== product type: SymmetricMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
// -> is identical with forwarding of Hermitian Matrix - GeneralMatrix products

//-- common interface ----------------------------------------------------------
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Side side,
       const ALPHA &alpha,
       const SymmetricMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, GeneralMatrix<MC> &C);

//-- symm
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Side side,
       const ALPHA &alpha, const SyMatrix<MA> &A, const GeMatrix<MB> &B,
       const BETA &beta, GeMatrix<MC> &C);


//== product type: TriangularMatrix - GeneralMatrix products

//-- forwarding ----------------------------------------------------------------
template <typename ALPHA, typename MA, typename MB>
    void
    mm(Side side, Transpose transA,
       const ALPHA &alpha, const MA &A, MB &&B);

//-- common interface ----------------------------------------------------------
template <typename ALPHA, typename MA, typename MB>
    void
    mm(Side side,
       Transpose transA, const ALPHA &alpha,
       const TriangularMatrix<MA> &A,
       GeneralMatrix<MB> &B);

//-- trmm
template <typename ALPHA, typename MA, typename MB>
    void
    mm(Side side,
       Transpose transA, const ALPHA &alpha, const TrMatrix<MA> &A,
       GeMatrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_MM_H
