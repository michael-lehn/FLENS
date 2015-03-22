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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_H
#define FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_H 1

#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//------------------------------------------------------------------------------
//
//  y = alpha*x
//
//  We evalute this with a scalSwitch
//  case 1: x is no closure
//  case 2: x is a closure

//
// case 1: x is vector but no closure
//
template <typename ALPHA, typename VX, typename VY>
    typename RestrictTo<!IsVectorClosure<Vector<VX>>::value,
             void>::Type
    copyScal(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y);

//
// case 2: x is a closure
//

//
// case 2 (a): x is a closure (conj)
//
template <typename ALPHA, typename VX, typename VY>
    typename RestrictTo<DefaultEval<VectorClosureOpConj<VX> >::value,
             void>::Type
    copyScal(const ALPHA &alpha, const VectorClosureOpConj<VX> &x,
             Vector<VY> &y);

//
// case 2 (b): x is a closure (general)
//
template <typename ALPHA, typename OP, typename VXL, typename VXR, typename VY>
    typename RestrictTo<VCDefaultEval<OP, VXL, VXR>::value,
             void>::Type
    copyScal(const ALPHA &alpha, const VectorClosure<OP, VXL, VXR> &x,
             Vector<VY> &y);

//------------------------------------------------------------------------------
//
//  B = alpha*op(A)
//
//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
// case 1: A is matrix but no closure
//
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<!IsMatrixClosure<Matrix<MA> >::value,
             void>::Type
    copyScal(Transpose trans, const ALPHA &alpha, const Matrix<MA> &A,
             Matrix<MB> &B);

//
// case 2: A is a closure
//

//
// case 2(a): A is a closure (trans)
//
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<DefaultEval<MatrixClosureOpTrans<MA> >::value,
             void>::Type
    copyScal(Transpose trans, const ALPHA &alpha,
             const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B);
//
// case 2(b): A is a closure (conj)
//
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<DefaultEval<MatrixClosureOpConj<MA> >::value,
             void>::Type
    copyScal(Transpose trans, const ALPHA &alpha,
             const MatrixClosureOpConj<MA> &A, Matrix<MB> &B);
//
// case 2(c): A is a closure (everything else)
//
template <typename ALPHA, typename OP, typename MAL, typename MAR, typename MB>
    typename RestrictTo<MCDefaultEval<OP, MAL, MAR>::value,
             void>::Type
    copyScal(Transpose trans, const ALPHA &alpha,
             const MatrixClosure<OP, MAL, MAR> &A, Matrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_H
