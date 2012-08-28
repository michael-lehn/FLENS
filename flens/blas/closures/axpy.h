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

#ifndef FLENS_BLAS_CLOSURES_AXPY_H
#define FLENS_BLAS_CLOSURES_AXPY_H 1

#include <cxxblas/typedefs.h>

#include <flens/blas/closures/closuretype.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- vector closures -----------------------------------------------------------
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

// y += scalar*x
template <typename ALPHA, typename T, typename VX, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpMult, ScalarValue<T>, VX> &x, Vector<VY> &y);

// y += x/scalar
template <typename ALPHA, typename VX, typename T, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<OpDiv, VX, ScalarValue<T> > &x, Vector<VY> &y);

// y += A*x
template <typename ALPHA, typename ML, typename VR, typename VY>
    typename RestrictTo<ClosureType<OpMult, ML, VR>::isMatrixVectorProduct,
             void>::Type
    axpy(const ALPHA &alpha,
         const VectorClosure<OpMult, ML, VR> &Ax, Vector<VY> &y);

// y += x*A
template <typename ALPHA, typename VL, typename MR, typename VY>
    typename RestrictTo<ClosureType<OpMult, VL, MR>::isVectorMatrixProduct,
             void>::Type
    axpy(const ALPHA &alpha,
         const VectorClosure<OpMult, VL, MR> &xA, Vector<VY> &y);

// y += <Unknown Closure>
template <typename ALPHA, typename Op, typename VL, typename VR, typename VY>
    void
    axpy(const ALPHA &alpha,
         const VectorClosure<Op, VL, VR> &x, Vector<VY> &y);

//-- matrix closures -----------------------------------------------------------
//
//  In the following comments op(X) denotes  x, X^T or X^H
//

//  B += alpha*op(A1 + A2)
template <typename ALPHA, typename ML, typename MR, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B);

//  B += alpha*op(A1 - A2)
template <typename ALPHA, typename ML, typename MR, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpSub, ML, MR> &A, Matrix<MB> &B);

//  B += scalar*op(A)
template <typename ALPHA, typename T, typename MA, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpMult, ScalarValue<T>, MA> &A, Matrix<MB> &B);

//  B += op(A)/scalar
template <typename ALPHA, typename MA, typename T, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpDiv, MA, ScalarValue<T> > &A, Matrix<MB> &B);

//  B += op(conjugate(A))
template <typename ALPHA, typename MA, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureOpConj<MA> &A, Matrix<MB> &B);

//  B += op(A^T)
template <typename ALPHA, typename MA, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B);

//  C += A*B
template <typename ALPHA, typename MA, typename MB, typename MC>
    typename RestrictTo<ClosureType<OpMult, MA, MB>::isMatrixMatrixProduct,
             void>::Type
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<OpMult,MA,MB> &AB, Matrix<MC> &C);

// B += <Unknown Closure>
template <typename ALPHA, typename Op, typename ML, typename MR, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<Op, ML, MR> &A, Matrix<MB> &B);

// B += Some Matrix
template <typename ALPHA, typename MA, typename MB>
    void
    axpy(Transpose trans, const ALPHA &alpha,
         const Matrix<MA> &A, Matrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_AXPY_H
