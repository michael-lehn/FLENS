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

#ifndef FLENS_BLAS_CLOSURES_DEBUGCLOSURE_H
#define FLENS_BLAS_CLOSURES_DEBUGCLOSURE_H 1

#include <flens/vectortypes/impl/densevector.h>
#include <flens/storage/fullstorage/hasfullstorage.h>

namespace flens {

//
//  If you implement a new matrix/vector type that supports views you might
//  need to specialize the following function
//

namespace DEBUGCLOSURE {

    template <typename X, typename Y>
        typename RestrictTo<!HasFullStorage<X>::value
                         || !HasFullStorage<Y>::value,
                 bool>::Type
        identical(const X &x, const Y &y);

    // TODO: Use a 'HasArrayStorage' instead of DenseVector
    template <typename VX, typename VY>
        bool
        identical(const DenseVector<VX> &x, const DenseVector<VY> &y);

    template <typename MA, typename MB>
        typename RestrictTo<HasFullStorage<MA>::value
                         && HasFullStorage<MB>::value,
                 bool>::Type
        identical(const MA &A, const MB &B);

} // namespace DEBUGCLOSURE


struct DebugClosure
{
    template <typename X, typename Y>
        static bool
        search(const X &x, const Y &y);

    template <typename Op, typename L, typename R, typename Y>
        static bool
        search(const VectorClosure<Op, L, R> &closure, const Y &y);

    template <typename Op, typename L, typename R, typename Y>
        static bool
        search(const MatrixClosure<Op, L, R> &closure, const Y &y);

};

} // namespace flens

#endif // FLENS_BLAS_CLOSURES_DEBUGCLOSURE_H
