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

#ifndef FLENS_STORAGE_ARRAY_ARRAYVIEW_TCC
#define FLENS_STORAGE_ARRAY_ARRAYVIEW_TCC 1

#include <cxxstd/cassert.h>
#include <cxxblas/level1/copy.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/array/array.h>
#include <flens/storage/array/arrayview.h>
#include <flens/storage/array/constarrayview.h>

namespace flens {

template <typename T, typename I, typename A>
ArrayView<T, I, A>::ArrayView(IndexType length, ElementType *data,
                              IndexType stride, IndexType firstIndex,
                              const Allocator &allocator)
    : _data(data),
      _allocator(allocator),
      _length(length),
      _stride(stride),
      _firstIndex(firstIndex)
{
    ASSERT(_length>=0);
    ASSERT(_stride>0);
}

template <typename T, typename I, typename A>
ArrayView<T, I, A>::ArrayView(const ArrayView &rhs)
    : _data(rhs._data),
      _allocator(rhs.allocator()),
      _length(rhs.length()),
      _stride(rhs.stride()),
      _firstIndex(rhs.firstIndex())
{
    ASSERT(_stride>0);
}

template <typename T, typename I, typename A>
template <typename RHS>
ArrayView<T, I, A>::ArrayView(RHS &rhs)
    : _data(rhs.data()),
      _allocator(rhs.allocator()),
      _length(rhs.length()),
      _stride(rhs.stride()),
      _firstIndex(rhs.firstIndex())
{
    ASSERT(_stride>0);
}

template <typename T, typename I, typename A>
ArrayView<T, I, A>::~ArrayView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I, typename A>
const typename ArrayView<T, I, A>::ElementType &
ArrayView<T, I, A>::operator()(IndexType index) const
{
#   ifndef NDEBUG
    if (lastIndex()>=firstIndex()) {
        ASSERT(index>=firstIndex());
        ASSERT(index<=lastIndex());
    }
#   endif

    return _data[_stride*(index-_firstIndex)];
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::ElementType &
ArrayView<T, I, A>::operator()(IndexType index)
{
#   ifndef NDEBUG
    if (lastIndex()>=firstIndex()) {
        ASSERT(index>=firstIndex());
        ASSERT(index<=lastIndex());
    }
#   endif

    return _data[_stride*(index-_firstIndex)];
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::IndexType
ArrayView<T, I, A>::firstIndex() const
{
    return _firstIndex;
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::IndexType
ArrayView<T, I, A>::lastIndex() const
{
    return _firstIndex+_length-1;
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::IndexType
ArrayView<T, I, A>::length() const
{
    return _length;
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::IndexType
ArrayView<T, I, A>::stride() const
{
    return _stride;
}

template <typename T, typename I, typename A>
const typename ArrayView<T, I, A>::ElementType *
ArrayView<T, I, A>::data() const
{
    return _data;
}

template <typename T, typename I, typename A>
typename ArrayView<T, I, A>::ElementType *
ArrayView<T, I, A>::data()
{
    return _data;
}

template <typename T, typename I, typename A>
const typename ArrayView<T, I, A>::Allocator &
ArrayView<T, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, typename I, typename A>
bool
ArrayView<T, I, A>::resize(IndexType DEBUG_VAR(length),
                           IndexType firstIndex,
                           const ElementType &)
{
    ASSERT(length==_length);

    changeIndexBase(firstIndex);
    return false;
}

template <typename T, typename I, typename A>
template <typename ARRAY>
bool
ArrayView<T, I, A>::resize(const ARRAY &rhs, const ElementType &value)
{
    return resize(rhs.length(), rhs.firstIndex(), value);
}

template <typename T, typename I, typename A>
bool
ArrayView<T, I, A>::fill(const ElementType &value)
{
    for (IndexType i=firstIndex(); i<=lastIndex(); ++i) {
        operator()(i) = value;
    }
    return true;
}

template <typename T, typename I, typename A>
void
ArrayView<T, I, A>::changeIndexBase(IndexType firstIndex)
{
    _firstIndex = firstIndex;
}

template <typename T, typename I, typename A>
const typename ArrayView<T, I, A>::ConstView
ArrayView<T, I, A>::view(IndexType from, IndexType to,
                         IndexType stride, IndexType firstViewIndex) const
{
    const IndexType length = (to-from)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ConstView(length,                // length
                         0,                     // data
                         stride*_stride,        // stride
                         firstViewIndex,        // firstIndex in view
                         allocator());          // allocator
    }
#   endif

    ASSERT(firstIndex()<=from);
    ASSERT(lastIndex()>=to);
    ASSERT(from<=to);
    ASSERT(stride>0);
    return ConstView(length,                // length
                     &operator()(from),     // data
                     stride*_stride,        // stride
                     firstViewIndex,        // firstIndex in view
                     allocator());          // allocator
}

template <typename T, typename I, typename A>
ArrayView<T, I, A>
ArrayView<T, I, A>::view(IndexType from, IndexType to,
                         IndexType stride, IndexType firstViewIndex)
{
    const IndexType length = (to-from)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ArrayView(length,                // length
                         0,                     // data
                         stride*_stride,        // stride
                         firstViewIndex,        // firstIndex in view
                         allocator());          // allocator
    } else {
        ASSERT(firstIndex()<=from);
        ASSERT(lastIndex()>=to);
        if (from>to) {
            std::cerr << "from = " << from << std::endl;
            std::cerr << "to =   " << to << std::endl;
            ASSERT(from<=to);
        }
    }
#   endif

    ASSERT(stride>0);
    return ArrayView(length,                // length
                     &operator()(from),     // data
                     stride*_stride,        // stride
                     firstViewIndex,        // firstIndex in view
                     allocator());          // allocator
}


//-- ArrayView specific functions ----------------------------------------------

//
//  fillRandom
//

template <typename T, typename I, typename A>
bool
fillRandom(ArrayView<T, I, A> &x)
{
    typedef typename ArrayView<T,I,A>::IndexType    IndexType;

    for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        x(i) = randomValue<T>();
    }
    return true;
}

} // namespace flens

#endif // FLENS_STORAGE_ARRAY_ARRAYVIEW_TCC
