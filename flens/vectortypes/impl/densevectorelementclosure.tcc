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

#ifndef FLENS_VECTORTYPES_IMPL_DENSEVECTORELEMENTCLOSURE_TCC
#define FLENS_VECTORTYPES_IMPL_DENSEVECTORELEMENTCLOSURE_TCC 1

namespace flens {

template <typename V>
DenseVectorElementClosure<V>::DenseVectorElementClosure(
                                                Vector &vector,
                                                IndexVariable &indexVariable)
    : _vector(vector), _indexVariable(indexVariable)
{
}

template <typename V>
void
DenseVectorElementClosure<V>::operator=(const ElementType &rhs)
{
    typename IndexVariable::ElementType &i = _indexVariable.value();

    for (i=_vector.firstIndex(); i<=_vector.lastIndex(); ++i) {
        value() = rhs;
    }
}

template <typename V>
template <typename S>
void
DenseVectorElementClosure<V>::operator=(const Scalar<S> &rhs)
{
    typename IndexVariable::ElementType &i = _indexVariable.value();

    for (i=_vector.firstIndex(); i<=_vector.lastIndex(); ++i) {
        value() = rhs.impl().value();
    }
}

template <typename V>
void
DenseVectorElementClosure<V>::operator=(const DenseVectorElementClosure &rhs)
{
    typename IndexVariable::ElementType &i = _indexVariable.value();

    for (i=_vector.firstIndex(); i<=_vector.lastIndex(); ++i) {
        value() = rhs.impl().value();
    }
}

template <typename V>
const typename DenseVectorElementClosure<V>::ElementType &
DenseVectorElementClosure<V>::value() const
{
    return _vector(_indexVariable.value());
}

template <typename V>
typename DenseVectorElementClosure<V>::ElementType &
DenseVectorElementClosure<V>::value()
{
    return _vector(_indexVariable.value());
}

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_DENSEVECTORELEMENTCLOSURE_TCC
