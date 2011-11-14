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

#ifndef FLENS_SCALAROPERATIONS_MINUS_TCC
#define FLENS_SCALAROPERATIONS_MINUS_TCC 1

namespace flens {

template <typename L, typename R>
const typename ScalarClosure<ScalarOpMinus, L, R>::ElementType
evalScalarClosure(const ScalarClosure<ScalarOpMinus, L, R> &exp)
{
    return exp.left().value() - exp.right().value();
}

//-- operator overloading
template <typename L, typename R>
const ScalarClosure<ScalarOpMinus,
                    typename L::Impl,
                    typename R::Impl>
operator-(const Scalar<L> &l, const Scalar<R> &r)
{
    typedef ScalarClosure<ScalarOpMinus, typename L::Impl, typename R::Impl> SC;
    return SC(l.impl(), r.impl());
}

template <typename ALPHA, typename R>
const typename 
RestrictTo<CompatibleScalar<ALPHA, R>::value, 
    ScalarClosure<ScalarOpMinus,
                  ScalarValue<ALPHA>,
                  typename R::Impl>
>::Type
operator-(const ALPHA &alpha, const Scalar<R> &r)
{
    typedef ScalarClosure<ScalarOpMinus,
                          ScalarValue<ALPHA>,
                          typename R::Impl>     SC;
    return SC(alpha, r.impl());
}

template <typename L, typename ALPHA>
const typename 
RestrictTo<CompatibleScalar<ALPHA, L>::value, 
    ScalarClosure<ScalarOpMinus,
                  typename L::Impl,
                  ScalarValue<ALPHA> >
>::Type
operator-(const Scalar<L> &l, const ALPHA &alpha)
{
    typedef ScalarClosure<ScalarOpMinus,
                          typename L::Impl,
                          ScalarValue<ALPHA> >  SC;
    return SC(l.impl(), alpha);
}

} // namespace flens

#endif // FLENS_SCALAROPERATIONS_MINUS_TCC
