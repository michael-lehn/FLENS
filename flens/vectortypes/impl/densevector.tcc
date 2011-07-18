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

#ifndef FLENS_VECTORTYPES_IMPL_DENSEVECTOR_TCC
#define FLENS_VECTORTYPES_IMPL_DENSEVECTOR_TCC 1

#include <flens/vectortypes/impl/dv/constelementclosure.tcc>
#include <flens/vectortypes/impl/dv/elementclosure.tcc>
#include <flens/vectortypes/impl/dv/initializer.tcc>

namespace flens {

template <typename A>
DenseVector<A>::DenseVector()
{
}

template <typename A>
DenseVector<A>::DenseVector(IndexType length)
    : _array(length)
{
}

template <typename A>
DenseVector<A>::DenseVector(IndexType length, IndexType firstIndex)
    : _array(length, firstIndex)
{
}

template <typename A>
DenseVector<A>::DenseVector(const Range<IndexType> &range)
    : _array(range.numTicks(), range.firstIndex())
{
}


template <typename A>
DenseVector<A>::DenseVector(const A &array)
    : _array(array)
{
}

template <typename A>
DenseVector<A>::DenseVector(const DenseVector &rhs)
    : Vector<DenseVector>(), _array(rhs._array)
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(const DenseVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(DenseVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(const Vector<RHS> &rhs)
{
    copy(rhs, *this);
}

// -- operators --------------------------------------------------------
template <typename A>
typename DenseVector<A>::Initializer
DenseVector<A>::operator=(const ElementType &value)
{
    engine().fill(value);
    return Initializer(*this, firstIndex());
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator=(const DenseVector &rhs)
{
    copy(rhs, *this);
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator=(const Vector<E> &rhs)
{
    copy(rhs, *this);
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator+=(const Vector<E> &rhs)
{
    axpy(ElementType(1), rhs, *this);
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator-=(const Vector<E> &rhs)
{
    axpy(ElementType(-1), rhs, *this);
    return *this;
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator+=(const ElementType &alpha)
{
	for (int i=firstIndex(); i<=lastIndex(); ++i) {
		(*this)(i) += alpha;
	}
    return *this;
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator-=(const ElementType &alpha)
{
	for (int i=firstIndex(); i<=lastIndex(); ++i) {
		(*this)(i) -= alpha;
	}
    return *this;
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator*=(const ElementType &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator/=(const ElementType &alpha)
{
    blas::scal(ElementType(1)/alpha, *this);
    return *this;
}

template <typename A>
const typename DenseVector<A>::ElementType &
DenseVector<A>::operator()(IndexType index) const
{
    return _array(index);
}

template <typename A>
typename DenseVector<A>::ElementType &
DenseVector<A>::operator()(IndexType index)
{
    return _array(index);
}

template <typename A>
template <typename S>
const densevector::ConstElementClosure<DenseVector<A>, typename Scalar<S>::Impl>
DenseVector<A>::operator()(const Scalar<S> &index) const
{
    typedef typename Scalar<S>::Impl ScalarImpl;
    typedef densevector::ConstElementClosure<DenseVector, ScalarImpl>  CEC;
    return CEC(*this, index.impl());
}

template <typename A>
const typename DenseVector<A>::ConstElementClosure
DenseVector<A>::operator()(const IndexVariable &index) const
{
    return ConstElementClosure(*this, index);
}

template <typename A>
typename DenseVector<A>::ElementClosure
DenseVector<A>::operator()(IndexVariable &index)
{
    return ElementClosure(*this, index);    
}

//-- views ---------------------------------------------------------------------

template <typename A>
typename DenseVector<A>::ConstView
DenseVector<A>::operator()(const Range<IndexType> &range) const
{
    ASSERT(range.firstIndex()>=firstIndex());
    ASSERT(range.lastIndex()<=lastIndex());
    return _array.view(range.firstIndex(), range.lastIndex(), range.stride(),
                       firstIndex());
}

template <typename A>
typename DenseVector<A>::View
DenseVector<A>::operator()(const Range<IndexType> &range)
{
    ASSERT(range.firstIndex()>=firstIndex());
    ASSERT(range.lastIndex()<=lastIndex());
    return _array.view(range.firstIndex(), range.lastIndex(), range.stride(),
                       firstIndex());
}

template <typename A>
typename DenseVector<A>::ConstView
DenseVector<A>::operator()(const Range<IndexType> &range,
                           IndexType _firstViewIndex) const
{
    ASSERT(range.firstIndex()>=firstIndex());
    ASSERT(range.lastIndex()<=lastIndex());
    return _array.view(range.firstIndex(), range.lastIndex(), range.stride(),
                       _firstViewIndex);
}

template <typename A>
typename DenseVector<A>::View
DenseVector<A>::operator()(const Range<IndexType> &range,
                           IndexType _firstViewIndex)
{
    ASSERT(range.firstIndex()>=firstIndex());
    ASSERT(range.lastIndex()<=lastIndex());
    return _array.view(range.firstIndex(), range.lastIndex(), range.stride(),
                       _firstViewIndex);
}

template <typename A>
typename DenseVector<A>::ConstView
DenseVector<A>::operator()(const Underscore<IndexType> &/*all*/,
                           IndexType _firstViewIndex) const
{
    return _array.view(_array.firstIndex(), _array.lastIndex(), _array.stride(),
                       _firstViewIndex);
}

template <typename A>
typename DenseVector<A>::View
DenseVector<A>::operator()(const Underscore<IndexType> &/*all*/,
                           IndexType _firstViewIndex)
{
    return _array.view(_array.firstIndex(), _array.lastIndex(), _array.stride(),
                       _firstViewIndex);
}

// -- methods ------------------------------------------------------------------

template <typename A>
Range<typename DenseVector<A>::IndexType>
DenseVector<A>::range() const
{
    return Range<IndexType>(_array.firstIndex(),_array.lastIndex());
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::firstIndex() const
{
    return _array.firstIndex();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::lastIndex() const
{
    return _array.lastIndex();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::length() const
{
    return _array.length();
}

template <typename A>
const A &
DenseVector<A>::engine() const
{
    return _array;
}

template <typename A>
A &
DenseVector<A>::engine()
{
    return _array;
}

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_DENSEVECTOR_TCC
