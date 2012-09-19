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

#ifndef FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGEVIEW_TCC
#define FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGEVIEW_TCC

#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/packedstorage/packedstorageview.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorageView<T, Order, I, A>::PackedStorageView(IndexType dim,
                                                     ElementType *data,
                                                     IndexType indexBase,
                                                     const Allocator &allocator)
    : _data(data),
      _allocator(allocator),
      _dim(dim),
      _indexBase(0)
{
    ASSERT(_dim>=0);

    changeIndexBase(indexBase);
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename ARRAY>
PackedStorageView<T, Order, I, A>::PackedStorageView(IndexType dim,
                                                     ARRAY &array,
                                                     IndexType indexBase,
                                                     const Allocator &allocator)
    : _data(array.data()),
      _allocator(allocator),
      _dim(dim),
      _indexBase(0)
{
    ASSERT((dim*(dim+1))/2<=array.length());
    ASSERT(_dim>=0);

    changeIndexBase(indexBase);
}

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorageView<T,Order,I,A>::PackedStorageView(const PackedStorageView &rhs)
    : _data(rhs._data),
      _allocator(rhs._allocator),
      _dim(rhs._dim),
      _indexBase(rhs._indexBase)
{
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
PackedStorageView<T, Order, I, A>::PackedStorageView(RHS &rhs)
    : _data(rhs.data()),
      _allocator(rhs.allocator()),
      _dim(rhs.dim()),
      _indexBase(0)
{
    changeIndexBase(rhs.indexBase());
}

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorageView<T, Order, I, A>::~PackedStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorageView<T, Order, I, A>::ElementType &
PackedStorageView<T, Order, I, A>::operator()(StorageUpLo  upLo,
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

    const IndexType i = row - _indexBase;
    const IndexType j = col - _indexBase;
    const IndexType n = _dim;

    if ((order==RowMajor) && (upLo==Upper)) {
        return _data[j+i*(2*n-i-1)/2];
    }
    if ((order==RowMajor) && (upLo==Lower)) {
        return _data[j+i*(i+1)/2];
    }
    if ((order==ColMajor) && (upLo==Upper)) {
        return _data[i+j*(j+1)/2];
    }
    return _data[i+j*(2*n-j-1)/2];
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorageView<T, Order, I, A>::ElementType &
PackedStorageView<T, Order, I, A>::operator()(StorageUpLo  upLo,
                                              IndexType    row,
                                              IndexType    col)
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

    const IndexType i = row - _indexBase;
    const IndexType j = col - _indexBase;
    const IndexType n = _dim;

    if ((order==RowMajor) && (upLo==Upper)) {
        return _data[j+i*(2*n-i-1)/2];
    }
    if ((order==RowMajor) && (upLo==Lower)) {
        return _data[j+i*(i+1)/2];
    }
    if ((order==ColMajor) && (upLo==Upper)) {
        return _data[i+j*(j+1)/2];
    }
    return _data[i+j*(2*n-j-1)/2];
}

//-- methods -------------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorageView<T, Order, I, A>::IndexType
PackedStorageView<T, Order, I, A>::indexBase() const
{
    return _indexBase;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorageView<T, Order, I, A>::IndexType
PackedStorageView<T, Order, I, A>::numNonZeros() const
{
    return _dim*(_dim+1)/IndexType(2);

}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorageView<T, Order, I, A>::IndexType
PackedStorageView<T, Order, I, A>::dim() const
{
    return _dim;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorageView<T, Order, I, A>::ElementType *
PackedStorageView<T, Order, I, A>::data() const
{
    return _data;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorageView<T, Order, I, A>::ElementType *
PackedStorageView<T, Order, I, A>::data()
{
    return _data;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorageView<T, Order, I, A>::Allocator &
PackedStorageView<T, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
PackedStorageView<T, Order, I, A>::resize(IndexType DEBUG_VAR(_dim),
                                          IndexType indexBase,
                                          const ElementType &)
{
    ASSERT(_dim==dim());
    changeIndexBase(indexBase);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
PackedStorageView<T, Order, I, A>::fill(const ElementType &value)
{
    std::fill_n(data(), numNonZeros(), value);
    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
PackedStorageView<T, Order, I, A>::changeIndexBase(IndexType indexBase)
{
    _indexBase = indexBase;
}

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGEVIEW_TCC
