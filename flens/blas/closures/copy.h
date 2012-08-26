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

#ifndef FLENS_BLAS_CLOSURES_COPY_H
#define FLENS_BLAS_CLOSURES_COPY_H 1

#include <cxxblas/typedefs.h>

#include <flens/blas/closures/closuretype.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//-- vector closures -----------------------------------------------------------
//  y = x1 + x2
template <typename VL, typename VR, typename VY>
    typename RestrictTo<!ClosureType<OpAdd, VL, VR>::isMatrixVectorProduct,
             void>::Type
    copy(const VectorClosure<OpAdd, VL, VR> &x, Vector<VY> &y);

//  y = x1 - x2
template <typename VL, typename VR, typename VY>
    typename RestrictTo<!ClosureType<OpSub, VL, VR>::isResidual,
             void>::Type
    copy(const VectorClosure<OpSub, VL, VR> &x, Vector<VY> &y);

//  y = alpha*x
template <typename T, typename VX, typename VY>
    void
    copy(const VectorClosure<OpMult, ScalarValue<T>, VX> &x, Vector<VY> &y);

//  y = x/alpha
template <typename VX, typename T, typename VY>
    void
    copy(const VectorClosure<OpDiv, VX, ScalarValue<T> > &x, Vector<VY> &y);

//  y = conjugate(x)
template <typename VX, typename VY>
    void
    copy(const VectorClosureOpConj<VX> &x, Vector<VY> &y);

// y = beta*y + alpha*op(A)*x
template <typename VL, typename MVR, typename VY>
    typename RestrictTo<ClosureType<OpAdd, VL, MVR>::isMatrixVectorProduct,
             void>::Type
    copy(const VectorClosure<OpAdd, VL, MVR> &ypAx, Vector<VY> &y);

// y = A*x
template <typename ML, typename VR, typename VY>
    typename RestrictTo<ClosureType<OpMult, ML, VR>::isMatrixVectorProduct,
             void>::Type
    copy(const VectorClosure<OpMult, ML, VR> &Ax, Vector<VY> &y);

// y = x*A
template <typename VL, typename MR, typename VY>
    typename RestrictTo<ClosureType<OpMult, VL, MR>::isVectorMatrixProduct,
             void>::Type
    copy(const VectorClosure<OpMult, VL, MR> &xA, Vector<VY> &y);

// y = b - A*x
template <typename L, typename R, typename VY>
    typename RestrictTo<ClosureType<OpSub, L, R>::isResidual,
             void>::Type
    copy(const VectorClosure<OpSub, L, R> &bmAx, Vector<VY> &y);

//  y = <Unknown Closure>
template <typename Op, typename VL, typename VR, typename VY>
    void
    copy(const VectorClosure<Op, VL, VR> &x, Vector<VY> &y);

//-- matrix closures -----------------------------------------------------------
//
//  In the following comments op(X) denotes  x, X^T or X^H
//

//  B = op(A1 + A2)
template <typename ML, typename MR, typename MB>
    typename RestrictTo<!ClosureType<OpAdd, ML, MR>::isMatrixMatrixProduct,
             void>::Type
    copy(Transpose trans, const MatrixClosure<OpAdd, ML, MR> &A, Matrix<MB> &B);

//  B = op(A1 - A2)
template <typename ML, typename MR, typename MB>
    void
    copy(Transpose trans, const MatrixClosure<OpSub, ML, MR> &A, Matrix<MB> &B);

//  B = alpha*op(A)
template <typename T, typename MA, typename MB>
    void
    copy(Transpose trans,
         const MatrixClosure<OpMult, ScalarValue<T>, MA> &A, Matrix<MB> &B);

//  B = op(A)/alpha
template <typename MA, typename T, typename MB>
    void
    copy(Transpose trans,
         const MatrixClosure<OpDiv, MA, ScalarValue<T> > &A, Matrix<MB> &B);

//  B = op(conjugate(A))
template <typename MA, typename MB>
    void
    copy(Transpose trans, const MatrixClosureOpConj<MA> &A, Matrix<MB> &B);

//  B = op(A^T)
template <typename MA, typename MB>
    void
    copy(Transpose trans, const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B);

//  C = beta*C + A*B
template <typename ML, typename MR, typename MC>
    typename RestrictTo<ClosureType<OpAdd, ML, MR>::isMatrixMatrixProduct,
             void>::Type
    copy(Transpose trans, const MatrixClosure<OpAdd, ML, MR> &CPAB,
         Matrix<MC> &C);

//  C = A*B
template <typename MA, typename MB, typename MC>
    typename RestrictTo<ClosureType<OpMult, MA, MB>::isMatrixMatrixProduct,
             void>::Type
    copy(Transpose trans, const MatrixClosure<OpMult,MA,MB> &AB, Matrix<MC> &C);

//  B = <Unknown Closure>
template <typename Op, typename ML, typename MR, typename MB>
    void
    copy(Transpose trans, const MatrixClosure<Op, ML, MR> &A, Matrix<MB> &B);

//-- hermitian matrices --------------------------------------------------------

template <typename MA, typename MB>
    void
    copy(Transpose trans, const HermitianMatrix<MA> &A, Matrix<MB> &B);

template <typename MA, typename MB>
    typename RestrictTo<!IsHermitianMatrix<MA>::value,
             void>::Type
    copy(Transpose trans, const Matrix<MA> &A, HermitianMatrix<MB> &B);

//-- symmetric matrices --------------------------------------------------------

template <typename MA, typename MB>
    void
    copy(Transpose trans, const SymmetricMatrix<MA> &A, Matrix<MB> &B);

template <typename MA, typename MB>
    typename RestrictTo<!IsSymmetricMatrix<MA>::value,
             void>::Type
    copy(Transpose trans, const Matrix<MA> &A, SymmetricMatrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_COPY_H
