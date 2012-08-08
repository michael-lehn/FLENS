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

#ifndef FLENS_BLAS_CLOSURES_DEBUGCLOSURE_TCC
#define FLENS_BLAS_CLOSURES_DEBUGCLOSURE_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/debugclosure.h>
#include <flens/vectortypes/vectortypes.h>

//
//  Default implementation and specializations of DEBUGCLOSURE::identical
//

namespace flens { namespace DEBUGCLOSURE {

template <typename X, typename Y>
typename RestrictTo<!HasFullStorage<X>::value || !HasFullStorage<Y>::value,
         bool>::Type
identical(const X &x, const Y &y)
{
    return ADDRESS(x)==ADDRESS(y);
}

//
// Two DenseVectors are identical if they reference the same memory
//
template <typename VX, typename VY>
bool
identical(const DenseVector<VX> &x, const DenseVector<VY> &y)
{
//
//  Quick return if possible
//
    if (ADDRESS(x)==ADDRESS(y)) {
        return true;
    }

    typedef typename DenseVector<VX>::ElementType  TX;
    typedef typename DenseVector<VY>::ElementType  TY;

    if (! IsSame<TX, TY>::value) {
        return false;
    }
//
//  Compare referenced memory
//
    if (RAWPOINTER(x.data())!=RAWPOINTER(y.data())) {
        return false;
    }
    if (x.stride()!=y.stride()) {
        return false;
    }
    if (x.length()!=y.length()) {
        return false;
    }
    return true;
}

//
// Two matrices with full storage are identical if they reference the
// same memory region
//
template <typename MA, typename MB>
typename RestrictTo<HasFullStorage<MA>::value && HasFullStorage<MB>::value,
         bool>::Type
identical(const MA &A, const MB &B)
{
//
//  Quick return if possible
//
    if (ADDRESS(A)==ADDRESS(B)) {
        return true;
    }

    typedef typename MA::ElementType  TA;
    typedef typename MB::ElementType  TB;

    if (! IsSame<TA, TB>::value) {
        return false;
    }

    if (MA::Engine::order!=MB::Engine::order) {
        return false;
    }
//
//  Compare referenced memory
//
    if (RAWPOINTER(A.data())!=RAWPOINTER(B.data())) {
        return false;
    }
    if (A.leadingDimension()!=B.leadingDimension()) {
        return false;
    }
    if (A.numRows()!=B.numRows()) {
        return false;
    }
    if (A.numCols()!=B.numCols()) {
        return false;
    }
    return true;
}

}} // namespace DEBUGCLOSURE, flens

//
//  Search mechanism for closures
//

namespace flens {

template <typename X, typename Y>
bool
DebugClosure::search(const X &x, const Y &y)
{
    return DEBUGCLOSURE::identical(x, y);
}

template <typename Op, typename L, typename R, typename Y>
bool
DebugClosure::search(const VectorClosure<Op, L, R> &closure, const Y &y)
{
    return search(closure.left(), y) || search(closure.right(), y);
}

template <typename Op, typename L, typename R, typename Y>
bool
DebugClosure::search(const MatrixClosure<Op, L, R> &closure, const Y &y)
{
    return search(closure.left(), y) || search(closure.right(), y);
}

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_DEBUGCLOSURE_TCC
