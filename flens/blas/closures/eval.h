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

#ifndef FLENS_BLAS_CLOSURES_EVAL_H
#define FLENS_BLAS_CLOSURES_EVAL_H 1

#include <cxxblas/typedefs.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

#ifdef USE_SPARSEBLAS_EXTENSIONS
#   include <extensions/flens/sparse_blas_flens.h>
#endif 

//-- blas entry points----------------------------------------------------------

namespace flens { 

//-- vector closures
template <typename VX, typename VY>
    void
    copy(const Vector<VX> &x, Vector<VY> &y);

template <typename ALPHA, typename VX, typename VY>
    void
    axpy(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y);

//-- matrix closures
template <typename MA, typename MB>
    void
    copy(cxxblas::Transpose trans,
         const Matrix<MA> &A, Matrix<MB> &B);

template <typename ALPHA, typename MA, typename MB>
    void
    axpy(cxxblas::Transpose trans,
         const ALPHA &alpha, const Matrix<MA> &A, Matrix<MB> &B);

} // namespace flens

namespace flens { namespace blas {

//-- vector closures -----------------------------------------------------------
//-- copy
// y = alpha*x
template <typename T, typename VX, typename VY>
    void
    copy(const VectorClosure<OpMult, ScalarValue<T>, VX> &ax, Vector<VY> &y);

// y = x1 + x2
template <typename VL, typename VR, typename VY>
    void
    copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y);

// y = x1 - x2
template <typename VL, typename VR, typename VY>
    void
    copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y);

// y = A*x
template <typename MA, typename VX, typename VY>
    void
    copy(const VectorClosure<OpMult, MA, VX> &Ax, Vector<VY> &y);

//-- axpy
// y += alpha*x
template <typename ALPHA, typename T, typename VX, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpMult, ScalarValue<T>, VX> &ax, Vector<VY> &y);

// y += x1 + x2
template <typename ALPHA, typename VL, typename VR, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y);

// y += x1 - x2
template <typename ALPHA, typename VL, typename VR, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y);

// y += A*x
template <typename ALPHA, typename MA, typename VX, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpMult, MA, VX> &Ax, Vector<VY> &y);

//-- mv
// y = x*A  ->  y = A'*x
template <typename ALPHA, typename VX, typename MA, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const Vector<VX> &x, const Matrix<MA> &A,
       const BETA &beta, Vector<VY> &y);

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const Matrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const HermitianMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const SymmetricMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

template <typename ALPHA, typename MA, typename VX, typename BETA, typename VY>
    void
    mv(cxxblas::Transpose trans,
       const ALPHA &alpha, const TriangularMatrix<MA> &A, const Vector<VX> &x,
       const BETA &beta, Vector<VY> &y);

//-- matrix closures -----------------------------------------------------------
//-- copy
// B = op(A)   (error-handle for unkown closure)
template <typename Op, typename MA1, typename MA2, typename MB>
    void
    copy(cxxblas::Transpose trans,
         const MatrixClosure<Op, MA1, MA2> &A, Matrix<MB> &B);

// B = A'
template <typename MA, typename MB>
    void
    copy(cxxblas::Transpose trans,
         const MatrixClosure<OpTrans, MA, MA> &At, Matrix<MB> &B);

// B = alpha*A
template <typename T, typename MA, typename MB>
    void
    copy(cxxblas::Transpose trans,
         const MatrixClosure<OpMult, ScalarValue<T>, MA> &aA, Matrix<MB> &B);

// B = A1 + A2
template <typename ML, typename MR, typename MB>
    void
    copy(cxxblas::Transpose trans,
         const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B);

//-- axpy
// B += op(A)   (error-handle for unkown closure)
template <typename ALPHA, typename Op, typename MA1, typename MA2, typename MB>
    void
    axpy(cxxblas::Transpose trans, const ALPHA &alpha,
         const MatrixClosure<Op, MA1, MA2> &A, Matrix<MB> &B);

// B += A'
template <typename ALPHA, typename MA, typename MB>
    void
    axpy(cxxblas::Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpTrans, MA, MA> &At, Matrix<MB> &B);

// B += alpha*A
template <typename ALPHA, typename T, typename MA, typename MB>
    void
    axpy(cxxblas::Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpMult, ScalarValue<T>, MA> &aA, Matrix<MB> &B);

// B += A1 + A2
template <typename ALPHA, typename ML, typename MR, typename MB>
    void
    axpy(cxxblas::Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_EVAL_H
