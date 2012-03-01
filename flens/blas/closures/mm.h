/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_MM_H
#define FLENS_BLAS_CLOSURES_MM_H 1

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//== GeneralMatrix - GeneralMatrix products ====================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const GeneralMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

//== TriangularMatrix - GeneralMatrix products =================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const TriangularMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const GeneralMatrix<MA> &A, const TriangularMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

//== SymmetricMatrix - GeneralMatrix products ==================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const SymmetricMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const GeneralMatrix<MA> &A, const SymmetricMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

//== HermitianMatrix - GeneralMatrix products ==================================
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const HermitianMatrix<MA> &A, const GeneralMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const GeneralMatrix<MA> &A, const HermitianMatrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

//== Matrix - Matrix products ==================================================
//
//  This gets called if everything else fails
//
#ifdef FLENS_DEBUG_CLOSURES

template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    mm(Transpose transA, Transpose transB, const ALPHA &alpha,
       const Matrix<MA> &A, const Matrix<MB> &B,
       const BETA &beta, Matrix<MC> &C);

#endif

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_MM_H

