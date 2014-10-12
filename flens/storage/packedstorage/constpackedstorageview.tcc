/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_TCC
#define FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_TCC 1

#include <flens/typedefs.h>
#include <flens/storage/packedstorage/constpackedstorageview.h>

namespace flens {

template <typename T, StorageOrder Order, typename I, typename A>
ConstPackedStorageView<T, Order, I, A>::ConstPackedStorageView(
                                                     IndexType dim,
                                                     const ElementType *data,
                                                     IndexType indexBase,
                                                     const Allocator &allocator)
    : data_(data), allocator_(allocator), dim_(dim),
      indexBase_(0)
{
    changeIndexBase(indexBase);
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename ARRAY>
ConstPackedStorageView<T, Order, I, A>::ConstPackedStorageView(
                                                     IndexType dim,
                                                     ARRAY &array,
                                                     IndexType indexBase,
                                                     const Allocator &allocator)
    : data_(array.data()),
      allocator_(allocator),
      dim_(dim),
      indexBase_(0)
{
    ASSERT(dim_>=0);
    ASSERT((dim_*(dim_+1))/2<=array.length());

    changeIndexBase(indexBase);
}

template <typename T, StorageOrder Order, typename I, typename A>
ConstPackedStorageView<T, Order, I, A>::ConstPackedStorageView(
                                              const ConstPackedStorageView &rhs)
    : data_(rhs.data_),
      allocator_(rhs.allocator_),
      dim_(rhs.dim_),
      indexBase_(rhs.indexBase_)
{
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
ConstPackedStorageView<T, Order, I, A>::ConstPackedStorageView(const RHS &rhs)
    : data_(rhs.data()),
      allocator_(rhs.allocator()),
      dim_(rhs.dim()),
      indexBase_(0)
{
    changeIndexBase(rhs.indexBase());
}

template <typename T, StorageOrder Order, typename I, typename A>
ConstPackedStorageView<T, Order, I, A>::~ConstPackedStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename ConstPackedStorageView<T, Order, I, A>::ElementType &
ConstPackedStorageView<T, Order, I, A>::operator()(StorageUpLo  upLo,
                                                   IndexType    row,
                                                   IndexType    col) const
{
#   ifndef NDEBUG
    ASSERT(row>=indexBase());
    ASSERT(col>=indexBase());
    ASSERT(row<=indexBase()+dim()-1);
    ASSERT(col<=indexBase()+dim()-1);

    if (dim()) {
        if (upLo==Lower) {
            ASSERT(row>=col);
        } else {
            ASSERT(col>=row);
        }
    }
#   endif

    const IndexType i = row - indexBase_;
    const IndexType j = col - indexBase_;
    const IndexType n = dim_;

    if ((order==RowMajor) && (upLo==Upper)) {
        return data_[j+i*(2*n-i-1)/2];
    }
    if ((order==RowMajor) && (upLo==Lower)) {
        return data_[j+i*(i+1)/2];
    }
    if ((order==ColMajor) && (upLo==Upper)) {
        return data_[i+j*(j+1)/2];
    }
    return data_[i+j*(2*n-j-1)/2];
}

//-- methods -------------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
typename ConstPackedStorageView<T, Order, I, A>::IndexType
ConstPackedStorageView<T, Order, I, A>::indexBase() const
{
    return indexBase_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename ConstPackedStorageView<T, Order, I, A>::IndexType
ConstPackedStorageView<T, Order, I, A>::numNonZeros() const
{
    return dim_*(dim_+1)/IndexType(2);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename ConstPackedStorageView<T, Order, I, A>::IndexType
ConstPackedStorageView<T, Order, I, A>::dim() const
{
    return dim_;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename ConstPackedStorageView<T, Order, I, A>::ElementType *
ConstPackedStorageView<T, Order, I, A>::data() const
{
    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename ConstPackedStorageView<T, Order, I, A>::Allocator &
ConstPackedStorageView<T, Order, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
ConstPackedStorageView<T, Order, I, A>::changeIndexBase(IndexType indexBase)
{
    indexBase_ = indexBase;

}

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_TCC
