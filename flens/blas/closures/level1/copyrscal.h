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

#ifndef FLENS_BLAS_CLOSURES_LEVEL1_COPYRSCAL_H
#define FLENS_BLAS_CLOSURES_LEVEL1_COPYRSCAL_H 1

#include <flens/blas/closures/tweaks/defaulteval.h>
#include <flens/blas/operators/operators.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/typedefs.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace blas {

//------------------------------------------------------------------------------
//
//  y = x/alpha
//
//  We evalute this with a scalSwitch
//  case 1: x is no closure
//  case 2: x is a closure

//
// case 1: x is vector but no closure
//
template <typename ALPHA, typename VX, typename VY>
    typename RestrictTo<IsVector<VX>::value
                    && !IsVectorClosure<VX>::value,
             void>::Type
    copyRScal(const ALPHA &alpha, const VX &x, Vector<VY> &y);

//
// case 2: x is a closure
//
template <typename ALPHA, typename VX, typename VY>
    typename RestrictTo<IsVectorClosure<VX>::value,
             void>::Type
    copyRScal(const ALPHA &alpha, const VX &x, Vector<VY> &y);

//------------------------------------------------------------------------------
//
//  B = op(A)/alpha
//
//  We evalute this with a scalSwitch
//  case 1: A is no closure
//  case 2: A is a closure

//
// case 1: A is matrix but no closure
//
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsMatrix<MA>::value
                    && !IsMatrixClosure<MA>::value,
             void>::Type
    copyRScal(Transpose trans, const ALPHA &alpha, const MA &A, Matrix<MB> &B);

//
// case 2: A is a closure
//
template <typename ALPHA, typename MA, typename MB>
    typename RestrictTo<IsMatrixClosure<MA>::value,
             void>::Type
    copyRScal(Transpose trans, const ALPHA &alpha, const MA &A, Matrix<MB> &B);

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_LEVEL1_COPYRSCAL_H
