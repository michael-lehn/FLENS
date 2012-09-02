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

#ifndef FLENS_BLAS_OPERATORS_OPMULT_H
#define FLENS_BLAS_OPERATORS_OPMULT_H 1

#include <flens/auxiliary/isconvertible.h>
#include <flens/auxiliary/restrictto.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/scalartypes/scalartypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens {

struct OpMult {};

//-- vector-vector products ----------------------------------------------------
// x'*y
template <typename VX, typename VY>
    typename Promotion<typename VX::Impl::ElementType,
                       typename VY::Impl::ElementType>::Type
    operator*(const Vector<VX> &x, const Vector<VY> &y);

//-- scalar-vector products ----------------------------------------------------
// alpha*x
template <typename ALPHA, typename V>
    const typename RestrictTo<!IsMatrix<ALPHA>::value
                           && !IsVector<ALPHA>::value,
          VectorClosure<OpMult, ScalarValue<ALPHA>, typename V::Impl> >::Type
    operator*(const ALPHA &alpha, const Vector<V> &x);

// x*alpha
template <typename ALPHA, typename V>
    const typename RestrictTo<!IsMatrix<ALPHA>::value
                           && !IsVector<ALPHA>::value,
          VectorClosure<OpMult, ScalarValue<ALPHA>, typename V::Impl> >::Type
    operator*(const Vector<V> &x, const ALPHA &alpha);

//-- scalar-matrix products ----------------------------------------------------
// alpha*A
template <typename ALPHA, typename M>
    const typename RestrictTo<!IsMatrix<ALPHA>::value
                           && !IsVector<ALPHA>::value,
          MatrixClosure<OpMult, ScalarValue<ALPHA>, typename M::Impl> >::Type
    operator*(const ALPHA &alpha, const Matrix<M> &A);

// A*alpha
template <typename ALPHA, typename M>
    const typename RestrictTo<!IsMatrix<ALPHA>::value
                           && !IsVector<ALPHA>::value,
          MatrixClosure<OpMult, ScalarValue<ALPHA>, typename M::Impl> >::Type
    operator*(const Matrix<M> &A, const ALPHA &alpha);

//-- matrix-vector products ----------------------------------------------------
// A*x
template <typename M, typename V>
    const VectorClosure<OpMult,
                        typename M::Impl,
                        typename V::Impl>
    operator*(const Matrix<M> &A, const Vector<V> &x);

// x*A
template <typename M, typename V>
    const VectorClosure<OpMult,
                        typename V::Impl,
                        typename M::Impl>
    operator*(const Vector<V> &x, const Matrix<M> &A);

//-- matrix-matrix products ----------------------------------------------------
// A*B
template <typename MA, typename MB>
    const MatrixClosure<OpMult,
                        typename MA::Impl,
                        typename MB::Impl>
    operator*(const Matrix<MA> &A, const Matrix<MB> &B);

} // namespace flens

#endif // FLENS_BLAS_OPERATORS_OPMULT_H
