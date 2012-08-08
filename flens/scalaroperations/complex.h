/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_SCALAROPERATIONS_COMPLEX_H
#define FLENS_SCALAROPERATIONS_COMPLEX_H 1

#include <cmath>
#include <flens/auxiliary/auxiliary.h>
#include <flens/scalartypes/scalartypes.h>

namespace flens {

struct ScalarOpComplex {};

//-- define the result of this operation is complex valued
template <typename L, typename R>
struct ElementType<ScalarClosure<ScalarOpComplex, L, R> >
{
    typedef typename L::ElementType                         LT;
    typedef typename R::ElementType                         RT;
    typedef std::complex<typename Promotion<LT, RT>::Type>  Type;
};

//-- eval: creates a complex value
template <typename L, typename R>
    const typename ScalarClosure<ScalarOpComplex, L, R>::ElementType
    evalScalarClosure(const ScalarClosure<ScalarOpComplex, L, R> &exp);

//-- function overloading
template <typename L, typename R>
    const ScalarClosure<ScalarOpComplex,
                        typename L::Impl,
                        typename R::Impl>
    Complex(const Scalar<L> &l, const Scalar<R> &r);

template <typename ALPHA, typename R>
    const typename
    RestrictTo<CompatibleScalar<ALPHA, R>::value,
        ScalarClosure<ScalarOpComplex,
                      ScalarValue<ALPHA>,
                      typename R::Impl>
    >::Type
    Complex(const ALPHA &alpha, const Scalar<R> &r);

template <typename ALPHA, typename L>
    const typename
    RestrictTo<CompatibleScalar<ALPHA, L>::value,
        ScalarClosure<ScalarOpComplex,
                      typename L::Impl,
                      ScalarValue<ALPHA> >
    >::Type
    Complex(const Scalar<L> &l, const ALPHA &alpha);

} // namespace flens

#endif // FLENS_SCALAROPERATIONS_COMPLEX_H
