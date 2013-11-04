/*
 *   Copyright (c) 2007, Michael Lehn
 *   Copyright (c) 2011, Michael Lehn
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

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/vectortypes/impl/densevector.h>
#include <flens/vectortypes/impl/dv/constelementclosure.tcc>
#include <flens/vectortypes/impl/dv/elementclosure.tcc>
#include <flens/vectortypes/impl/dv/initializer.tcc>

namespace flens {

template <typename A>
DenseVector<A>::DenseVector()
    : _stride(_array.stride())
{
}

template <typename A>
DenseVector<A>::DenseVector(IndexType length)
    : _array(length), _stride(_array.stride())
{
    ASSERT(length>=0);
}

template <typename A>
DenseVector<A>::DenseVector(IndexType length, IndexType firstIndex)
    : _array(length, firstIndex), _stride(_array.stride())
{
    ASSERT(length>=0);
}

template <typename A>
DenseVector<A>::DenseVector(const Range<IndexType> &range)
    : _array(range.numTicks(), range.firstIndex()), _stride(_array.stride())
{
    ASSERT(range.stride()>0);
}

template <typename A>
DenseVector<A>::DenseVector(const Engine &engine, bool reverse)
    : _array(engine), _stride(reverse ? -_array.stride() : _array.stride())
{
}

template <typename A>
DenseVector<A>::DenseVector(const DenseVector &rhs)
    : Vector<DenseVector>(), _array(rhs._array), _stride(rhs._stride)
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(const DenseVector<RHS> &rhs)
    : _array(rhs.engine()), _stride(_array.stride())
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(DenseVector<RHS> &rhs)
    : _array(rhs.engine()), _stride(_array.stride())
{
}

template <typename A>
template <typename RHS>
DenseVector<A>::DenseVector(const Vector<RHS> &rhs)
    : _stride(_array.stride())
{
    assign(rhs, *this);
}

//-- operators -----------------------------------------------------------------
template <typename A>
typename DenseVector<A>::Initializer
DenseVector<A>::operator=(const ElementType &value)
{
    _array.fill(value);
    return Initializer(*this, firstIndex());
}

template <typename A>
DenseVector<A> &
DenseVector<A>::operator=(const DenseVector &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator=(const Vector<E> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator+=(const Vector<E> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename A>
template <typename E>
DenseVector<A> &
DenseVector<A>::operator-=(const Vector<E> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename A>
template <typename T>
typename RestrictTo<IsConvertible<T,typename A::ElementType>::value,
                    DenseVector<A> >::Type &
DenseVector<A>::operator+=(const T &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) += alpha;
    }
    return *this;
}

template <typename A>
template <typename T>
typename RestrictTo<IsConvertible<T,typename A::ElementType>::value,
                    DenseVector<A> >::Type &
DenseVector<A>::operator-=(const T &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) -= alpha;
    }
    return *this;
}

template <typename A>
template <typename T>
typename RestrictTo<IsConvertible<T,typename A::ElementType>::value,
                    DenseVector<A> >::Type &
DenseVector<A>::operator*=(const T &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename A>
template <typename T>
typename RestrictTo<IsConvertible<T,typename A::ElementType>::value,
                    DenseVector<A> >::Type &
DenseVector<A>::operator/=(const T &alpha)
{
    blas::rscal(alpha, *this);
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
const typename DenseVector<A>::ConstView
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
const typename DenseVector<A>::ConstView
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
const typename DenseVector<A>::ConstView
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

template <typename A>
const typename DenseVector<A>::ConstView
DenseVector<A>::reverse() const
{
    return ConstView(_array, !reversed());
}

template <typename A>
typename DenseVector<A>::View
DenseVector<A>::reverse()
{
    return View(_array, !reversed());
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
    if (_stride<0) {
        return _array.lastIndex();
    }
    return _array.firstIndex();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::lastIndex() const
{
    if (_stride<0) {
        return _array.firstIndex();
    }
    return _array.lastIndex();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::length() const
{
    return _array.length();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::inc() const
{
    if (_stride<0) {
        return -1;
    }
    return 1;
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::endIndex() const
{
    return lastIndex() + inc();
}

template <typename A>
const typename DenseVector<A>::ElementType *
DenseVector<A>::data() const
{
    return _array.data();
}

template <typename A>
typename DenseVector<A>::ElementType *
DenseVector<A>::data()
{
    return _array.data();
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::stride() const
{
#   ifndef NDEBUG
    if (!reversed()) {
        ASSERT(_stride==_array.stride());
    } else {
        ASSERT(_stride==-_array.stride());
    }
#   endif

    return _stride;
}

template <typename A>
typename DenseVector<A>::IndexType
DenseVector<A>::indexBase() const
{
    return _array.firstIndex();
}

template <typename A>
template <typename RHS>
bool
DenseVector<A>::resize(const DenseVector<RHS> &rhs, const ElementType &value)
{
    bool resized = _array.resize(rhs.engine(), value);
    return resized;
}

template <typename A>
bool
DenseVector<A>::resize(IndexType length, IndexType firstIndex,
                       const ElementType &value)
{
    return _array.resize(length, firstIndex, value);
}

template <typename A>
bool
DenseVector<A>::fill(const ElementType &value)
{
    return _array.fill(value);
}

template <typename A>
void
DenseVector<A>::changeIndexBase(IndexType firstIndex)
{
    _array.changeIndexBase(firstIndex);
}

// -- implementation -----------------------------------------------------------
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

template <typename A>
bool
DenseVector<A>::reversed() const
{
    return (_stride<0);
}

//-- DenseVector specific functions --------------------------------------------

//
//  real
//

template <typename VZ>
typename RealVector<DenseVector<VZ> >::ConstView
real(const DenseVector<VZ> &z)
{
    typedef typename RealVector<DenseVector<VZ> >::ConstViewEngine  Engine;
    typedef typename VZ::ElementType                                T;
    typedef typename Engine::ElementType                            PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<const PT *>(z.data()),
                  z.stride()*sizeof(T)/sizeof(PT),
                  z.firstIndex());
}

template <typename VZ>
typename RestrictTo<IsDenseVector<VZ>::value,
         typename RealVector<typename RemoveRef<VZ>::Type>::View>::Type
real(VZ &&z)
{
    typedef typename RemoveRef<VZ>::Type                VectorZ;
    typedef typename RealVector<VectorZ>::ViewEngine    Engine;
    typedef typename VectorZ::ElementType               T;
    typedef typename Engine::ElementType                PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<PT *>(z.data()),
                  z.stride()*sizeof(T)/sizeof(PT),
                  z.firstIndex());
}


//
//  imag
//

template <typename VZ>
typename ImagVector<DenseVector<VZ> >::ConstView
imag(const DenseVector<VZ> &z)
{
    typedef typename ImagVector<DenseVector<VZ> >::ConstViewEngine  Engine;
    typedef typename VZ::ElementType                                T;
    typedef typename Engine::ElementType                            PT;

    static_assert(sizeof(T)==2*sizeof(PT), "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<const PT *>(z.data()) + 1,
                  2*z.stride(),
                  z.firstIndex());
}

template <typename VZ>
typename RestrictTo<IsDenseVector<VZ>::value,
         typename ImagVector<typename RemoveRef<VZ>::Type>::View>::Type
imag(VZ &&z)
{
    typedef typename RemoveRef<VZ>::Type                VectorZ;
    typedef typename ImagVector<VectorZ>::ViewEngine    Engine;
    typedef typename VectorZ::ElementType               T;
    typedef typename Engine::ElementType                PT;

    static_assert(sizeof(T)%sizeof(PT)==0, "Are you kidding me?");

    return Engine(z.length(),
                  reinterpret_cast<PT *>(z.data()) + 1,
                  2*z.stride(),
                  z.firstIndex());
}


//
//  fillRandom
//

template <typename VX>
typename RestrictTo<IsDenseVector<VX>::value,
         bool>::Type
fillRandom(VX &&x)
{
    return fillRandom(x.engine());
}


} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_DENSEVECTOR_TCC
