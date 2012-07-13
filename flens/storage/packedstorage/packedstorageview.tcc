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

// #include <flens/storage/fullstorage/trapezoidalfill.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
PackedStorageView<T, UpLo, Order, I, A>::PackedStorageView(IndexType dim,
                                                          ElementType *data,
                                                          IndexType firstIndex,
                                                          const Allocator &allocator)
    : _data(data),
      _allocator(allocator),
      _dim(dim),
      _firstIndex(firstIndex)
{
    ASSERT(_dim>=0);

    changeIndexBase(firstIndex);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
template <typename ARRAY>
PackedStorageView<T, UpLo, Order, I, A>::PackedStorageView(IndexType dim,
                                                           ARRAY &array,
                                                           IndexType firstIndex,
                                                           const Allocator &allocator)
    : _data(array.data()),
      _allocator(allocator),
      _dim(dim),
      _firstIndex(0)
{
    ASSERT((dim*(dim+1))/2<=array.length());
    ASSERT(_dim>=0);

    changeIndexBase(firstIndex);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
PackedStorageView<T, UpLo, Order, I, A>::PackedStorageView(const PackedStorageView &rhs)
    : _data(rhs._data),
      _allocator(rhs._allocator),
      _dim(rhs._dim),
      _firstIndex(rhs._firstIndex)
{
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
template <typename RHS>
PackedStorageView<T, UpLo, Order, I, A>::PackedStorageView(RHS &rhs)
    : _data(rhs.data()),
      _allocator(rhs.allocator()),
      _dim(rhs.dim()),
      _firstIndex(0)
{
    changeIndexBase(rhs.firstIndex());
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
PackedStorageView<T, UpLo, Order, I, A>::~PackedStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
const typename PackedStorageView<T, UpLo, Order, I, A>::ElementType &
PackedStorageView<T, UpLo, Order, I, A>::operator()(IndexType row, IndexType col) const
{
    assert(row>=firstIndex());
    assert(col>=firstIndex());
    assert(row<=lastIndex());
    assert(col<=lastIndex());

    if (UpLo==Lower) {
        assert(row>=col);
    } else {
        assert(col>=row);
    }

    const IndexType i = row - _firstIndex;
    const IndexType j = col - _firstIndex;
    const IndexType n = _dim;

    if ((order==RowMajor) && (UpLo==Upper)) {
        return _data[j+i*(2*n-i-1)/2];
    }
    if ((order==RowMajor) && (UpLo==Lower)) {
        return _data[j+i*(i+1)/2];
    }
    if ((order==ColMajor) && (UpLo==Upper)) {
        return _data[i+j*(j+1)/2];
    }
    return _data[i+j*(2*n-j-1)/2];
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
typename PackedStorageView<T, UpLo, Order, I, A>::ElementType &
PackedStorageView<T, UpLo, Order, I, A>::operator()(IndexType row, IndexType col)
{
    assert(row>=firstIndex());
    assert(col>=firstIndex());
    assert(row<=lastIndex());
    assert(col<=lastIndex());

    if (UpLo==Lower) {
        assert(row>=col);
    } else {
        assert(col>=row);
    }

    const IndexType i = row - _firstIndex;
    const IndexType j = col - _firstIndex;
    const IndexType n = _dim;

    if ((order==RowMajor) && (UpLo==Upper)) {
        return _data[j+i*(2*n-i-1)/2];
    }
    if ((order==RowMajor) && (UpLo==Lower)) {
        return _data[j+i*(i+1)/2];
    }
    if ((order==ColMajor) && (UpLo==Upper)) {
        return _data[i+j*(j+1)/2];
    }
    return _data[i+j*(2*n-j-1)/2];
}

//-- methods -------------------------------------------------------------------

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
typename PackedStorageView<T, UpLo, Order, I, A>::IndexType
PackedStorageView<T, UpLo, Order, I, A>::firstIndex() const
{
    return _firstIndex;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
typename PackedStorageView<T, UpLo, Order, I, A>::IndexType
PackedStorageView<T, UpLo, Order, I, A>::lastIndex() const
{
    return _firstIndex+_dim-1;

}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
typename PackedStorageView<T, UpLo, Order, I, A>::IndexType
PackedStorageView<T, UpLo, Order, I, A>::dim() const
{
    return _dim;
}


template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
const typename PackedStorageView<T, UpLo, Order, I, A>::ElementType *
PackedStorageView<T, UpLo, Order, I, A>::data() const
{
#   ifndef NDEBUG
    if (dim()==0) {
        return 0;
    }
#   endif
    return &(this->operator()(_firstIndex, _firstIndex));
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
typename PackedStorageView<T, UpLo, Order, I, A>::ElementType *
PackedStorageView<T, UpLo, Order, I, A>::data()
{
#   ifndef NDEBUG
    if (dim()==0) {
        return 0;
    }
#   endif

    return &(this->operator()(_firstIndex, _firstIndex));
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
const typename PackedStorageView<T, UpLo, Order, I, A>::Allocator &
PackedStorageView<T, UpLo, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
bool
PackedStorageView<T, UpLo, Order, I, A>::resize(IndexType DEBUG_VAR(_dim),
                                        IndexType firstIndex,
                                        const ElementType &)
{
    ASSERT(_dim==dim());
    changeIndexBase(firstIndex);
    return false;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
template <typename FS>
bool
PackedStorageView<T, UpLo, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.firstIndex(), rhs.dim(),
                  value);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
bool
PackedStorageView<T, UpLo, Order, I, A>::fill(const ElementType &value)
{
    std::fill_n(data(), 0.5*(dim()+1)*dim(), value);
    return true;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order,  typename I, typename A>
void
PackedStorageView<T, UpLo, Order, I, A>::changeIndexBase(IndexType firstIndex)
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (dim()==0)  {
        _firstIndex = firstIndex;
        return;
    }
#   endif

    _firstIndex = firstIndex;
}




} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGEVIEW_TCC
