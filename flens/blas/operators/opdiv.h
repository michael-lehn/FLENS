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

#ifndef FLENS_BLAS_OPERATORS_OPDIV_H
#define FLENS_BLAS_OPERATORS_OPDIV_H 1

#include <flens/aux/isconvertible.h>
#include <flens/aux/restrictto.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/scalartypes/scalartypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens {

struct OpDiv {};

//
//  x/alpha
//
//  The declaration of this function certainly looks complicated.  But what
//  the RestrictTo trait does is ensuring that alpha is a scalar and
//   - either alpha can be converted to the element type of x
//   - or the element type of x can beconverted to alpha.
//  If this is the case the return type is
//
//      const VectorClosure<OpDiv, typename V::Impl, ScalarValue<ALPHA> >
//
template <typename ALPHA, typename V>
    const typename RestrictTo<
            IsConvertible<ALPHA, typename V::Impl::ElementType>::value
         || IsConvertible<typename V::Impl::ElementType, ALPHA>::value,
            VectorClosure<OpDiv, typename V::Impl, ScalarValue<ALPHA> >
        >::Type
    operator/(const Vector<V> &x, const ALPHA &alpha);

//
//  A/alpha
//
template <typename ALPHA, typename M>
    const typename RestrictTo<
            IsConvertible<ALPHA, typename M::Impl::ElementType>::value
         || IsConvertible<typename M::Impl::ElementType, ALPHA>::value,
            MatrixClosure<OpDiv, typename M::Impl, ScalarValue<ALPHA> >
        >::Type
    operator/(const Matrix<M> &x, const ALPHA &alpha);

} // namespace flens

#endif // FLENS_BLAS_OPERATORS_OPDIV_H
