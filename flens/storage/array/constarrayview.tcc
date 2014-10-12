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

#ifndef FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_TCC
#define FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_TCC 1

#include <cxxstd/cassert.h>
#include <flens/storage/array/array.h>
#include <flens/storage/array/arrayview.h>
#include <flens/storage/array/constarrayview.h>


namespace flens {

template <typename T, typename I, typename A>
ConstArrayView<T, I, A>::ConstArrayView(IndexType length,
                                        const ElementType *data,
                                        IndexType stride,
                                        IndexType firstIndex,
                                        const Allocator &allocator)
    : data_(data),
      allocator_(allocator),
      length_(length),
      stride_(stride),
      firstIndex_(firstIndex)
{
    ASSERT(length_>=0);
    ASSERT(stride_>0);
}

template <typename T, typename I, typename A>
ConstArrayView<T, I, A>::ConstArrayView(const ConstArrayView &rhs)
    : data_(rhs.data_),
      allocator_(rhs.allocator()),
      length_(rhs.length_),
      stride_(rhs.stride_),
      firstIndex_(rhs.firstIndex_)
{
    ASSERT(stride_>0);
}

template <typename T, typename I, typename A>
template <typename RHS>
ConstArrayView<T, I, A>::ConstArrayView(const RHS &rhs)
    : data_(rhs.data()),
      allocator_(rhs.allocator()),
      length_(rhs.length()),
      stride_(rhs.stride()),
      firstIndex_(rhs.firstIndex())
{
    ASSERT(stride_>0);
}

template <typename T, typename I, typename A>
ConstArrayView<T, I, A>::~ConstArrayView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I, typename A>
const typename ConstArrayView<T, I, A>::ElementType &
ConstArrayView<T, I, A>::operator()(IndexType index) const
{
    ASSERT(index>=firstIndex());
    ASSERT(index<=lastIndex());
    return data_[stride_*(index-firstIndex_)];
}

//-- methods -------------------------------------------------------------------

template <typename T, typename I, typename A>
typename ConstArrayView<T, I, A>::IndexType
ConstArrayView<T, I, A>::firstIndex() const
{
    return firstIndex_;
}

template <typename T, typename I, typename A>
typename ConstArrayView<T, I, A>::IndexType
ConstArrayView<T, I, A>::lastIndex() const
{
    return firstIndex_+length_-1;
}

template <typename T, typename I, typename A>
typename ConstArrayView<T, I, A>::IndexType
ConstArrayView<T, I, A>::length() const
{
    return length_;
}

template <typename T, typename I, typename A>
typename ConstArrayView<T, I, A>::IndexType
ConstArrayView<T, I, A>::stride() const
{
    return stride_;
}

template <typename T, typename I, typename A>
const typename ConstArrayView<T, I, A>::ElementType *
ConstArrayView<T, I, A>::data() const
{
    return data_;
}

template <typename T, typename I, typename A>
const typename ConstArrayView<T, I, A>::Allocator &
ConstArrayView<T, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, typename I, typename A>
const ConstArrayView<T, I, A>
ConstArrayView<T, I, A>::view(IndexType from, IndexType to,
                              IndexType stride, IndexType firstViewIndex) const
{
    const IndexType length = (to-from)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ConstArrayView<T, I, A>(length,          // length
                                       0,               // data
                                       stride*stride_,  // stride
                                       firstViewIndex,  // firstIndex in view
                                       allocator());    // allocator
    }
#   endif

    ASSERT(firstIndex()<=from);
    ASSERT(lastIndex()>=to);
    ASSERT(from<=to);
    ASSERT(stride>=1);
    return ConstArrayView<T, I, A>(length,              // length
                                   &operator()(from),   // data
                                   stride*stride_,      // stride
                                   firstViewIndex,      // firstIndex in view
                                   allocator());        // allocator
}

template <typename T, typename I, typename A>
void
ConstArrayView<T, I, A>::changeIndexBase(IndexType firstIndex)
{
    firstIndex_ = firstIndex;
}

} // namespace flens

#endif // FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_TCC
