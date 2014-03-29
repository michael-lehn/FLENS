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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_TCC
#define FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

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
copyScal(const ALPHA &alpha, const Vector<VX> &x, Vector<VY> &y)
{
    using namespace DEBUGCLOSURE;

    if (identical(x.impl(), y.impl())) {
        scal(alpha, y.impl());
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        const typename VX::Impl::ElementType  Zero(0);
        axpby(alpha, x.impl(), Zero, y.impl());
    }
}

template <typename ALPHA, typename VX, typename VY>
typename RestrictTo<DefaultEval<VectorClosureOpConj<VX> >::value,
         void>::Type
copyScal(const ALPHA &alpha, const VectorClosureOpConj<VX> &x, Vector<VY> &y)
{
    using namespace DEBUGCLOSURE;
    if (identical(x.left(), y.impl())) {
        copyConj(x.left(), y.impl());
        scal(alpha, y.impl());
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        const typename VX::Impl::ElementType  Zero(0);
        acxpby(alpha, x.left(), Zero, y.impl());
    }
}
//
// case 2 (b): x is a closure (general)
//
template <typename ALPHA, typename OP, typename VXL, typename VXR, typename VY>
typename RestrictTo<VCDefaultEval<OP, VXL, VXR>::value,
         void>::Type
copyScal(const ALPHA &alpha, const VectorClosure<OP, VXL, VXR> &x,
         Vector<VY> &y)
{
    copy(x, y.impl());
    scal(alpha, y.impl());
}

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
         Matrix<MB> &B)
{
    using namespace DEBUGCLOSURE;

    if (identical(A.impl(), B.impl())) {
        scal(alpha, B.impl());
        if (trans!=NoTrans) {
            //cotr(trans, B.impl());
        }
    } else {
//
//      Zero should have the same type as elements of y so that CBLAS
//      gets called if possible.
//
        const typename Matrix<MA>::Impl::ElementType  Zero(0);
        axpby(trans, alpha, A.impl(), Zero, B.impl());
        //scal(Zero, B.impl());
        //axpy(trans, alpha, A.impl(), B.impl());
    }
}

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
         const MatrixClosureOpTrans<MA> &A, Matrix<MB> &B)
{
    trans = Transpose(trans^Trans);
    copyScal(trans, alpha, A.left(), B.impl());
}

//
// case 2(b): A is a closure (conj)
//
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<DefaultEval<MatrixClosureOpConj<MA> >::value,
         void>::Type
copyScal(Transpose trans, const ALPHA &alpha,
         const MatrixClosureOpConj<MA> &A, Matrix<MB> &B)
{
    trans = Transpose(trans^Conj);
    copyScal(trans, alpha, A.left(), B.impl());
}

//
// case 2(c): A is a closure (everything else)
//
template <typename ALPHA, typename Op, typename MAL, typename MAR, typename MB>
typename RestrictTo<MCDefaultEval<Op, MAL, MAR>::value,
         void>::Type
copyScal(Transpose trans, const ALPHA &alpha,
         const MatrixClosure<Op, MAL, MAR> &A, Matrix<MB> &B)
{
    copy(trans, A, B.impl());
    scal(alpha, B.impl());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPYSCAL_TCC
