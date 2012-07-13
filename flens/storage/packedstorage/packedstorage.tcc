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

#ifndef FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_TCC
#define FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_TCC

#include <cxxblas/level1/copy.h>
#include <flens/typedefs.h>

namespace flens {

//= Constructors

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
PackedStorage<T, UpLo, Order, I, A>::PackedStorage()
    :  _data(0),
       _dim(0),
       _firstIndex(I::defaultIndexBase)
{
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
PackedStorage<T, UpLo, Order, I, A>::PackedStorage(IndexType dim,
                                                   IndexType firstIndex,
                                                   const ElementType &value,
                                                   const Allocator &allocator)
    : _data(0), _allocator(allocator),
      _dim(dim),
      _firstIndex(firstIndex)
{
    ASSERT(_dim>=0);

    _allocate(value);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
PackedStorage<T, UpLo, Order, I, A>::PackedStorage(const PackedStorage &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _dim(rhs.dim()),
      _firstIndex(rhs.firstIndex())
{
    _allocate(ElementType());

    cxxblas::copy(rhs.dim(), rhs.data(), 1, data(), 1);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
template <typename RHS>
PackedStorage<T, UpLo, Order, I, A>::PackedStorage(const RHS &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _dim(rhs.dim()),
      _firstIndex(rhs.firstIndex())
{
    _allocate(ElementType());

    cxxblas::copy(rhs.dim(), rhs.data(), 1, data(), 1);

}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
PackedStorage<T, UpLo, Order, I, A>::~PackedStorage()
{
    _release();
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, UpLo, Order, I, A>::ElementType &
PackedStorage<T, UpLo, Order, I, A>::operator()(IndexType row, IndexType col) const
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

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, UpLo, Order, I, A>::ElementType &
PackedStorage<T, UpLo, Order, I, A>::operator()(IndexType row, IndexType col)
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

//-- Methods -------------------------------------------------------------------
template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, UpLo, Order, I, A>::IndexType
PackedStorage<T, UpLo, Order, I, A>::firstIndex() const
{
    return _firstIndex;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, UpLo, Order, I, A>::IndexType
PackedStorage<T, UpLo, Order, I, A>::lastIndex() const
{
    return _firstIndex+_dim-1;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, UpLo, Order, I, A>::IndexType
PackedStorage<T, UpLo, Order, I, A>::dim() const
{
    return _dim;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, UpLo, Order, I, A>::ElementType *
PackedStorage<T, UpLo, Order, I, A>::data() const
{
#   ifndef NDEBUG
    if (dim()==0) {
        return 0;
    }
#   endif

    ASSERT(_data);
    return &(this->operator()(_firstIndex, _firstIndex));
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, UpLo, Order, I, A>::ElementType *
PackedStorage<T, UpLo, Order, I, A>::data()
{
#   ifndef NDEBUG
    if (dim()==0)  {
        return 0;
    }
#   endif

    ASSERT(_data);
    return &(this->operator()(_firstIndex, _firstIndex));
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, UpLo, Order, I, A>::Allocator &
PackedStorage<T, UpLo, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
bool
PackedStorage<T, UpLo, Order, I, A>::resize(IndexType dim,
                                    IndexType firstIndex,
                                    const ElementType &value)
{
    if (_dim!=dim) {
        _release();
        _dim = dim;
        _firstIndex = firstIndex;
        _allocate(value);
        return true;
    }
    changeIndexBase(firstIndex);
    return false;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
bool
PackedStorage<T, UpLo, Order, I, A>::resize(const Range<IndexType> &rac,
                                            const ElementType &value)
{
    if (_dim!=rac.length()){
        _release();
        _dim = rac.length();
        _firstIndex = rac.firstIndex();
        _allocate(value);
        return true;
    }
    changeIndexBase(rac.firstIndex());
    return false;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
template <typename FS>
bool
PackedStorage<T, UpLo, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.dim(),
                  rhs.firstIndex(),
                  value);
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
bool
PackedStorage<T, UpLo, Order, I, A>::fill(const ElementType &value)
{
    ASSERT(_data);
    std::fill_n(data(), 0.5*(dim()+1)*dim(), value);
    return true;
}


template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, UpLo, Order, I, A>::changeIndexBase(IndexType firstIndex)
{
    _firstIndex = firstIndex;
}




//-- Private Methods -----------------------------------------------------------

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, UpLo, Order, I, A>::_setIndexBase(IndexType firstIndex)
{
    _firstIndex = firstIndex;
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, UpLo, Order, I, A>::_raw_allocate()
{
    ASSERT(!_data);
    ASSERT(_dim>0);

    _data = _allocator.allocate(0.5*(_dim+1)*_dim);
#ifndef NDEBUG
    ElementType *p = _data;
#endif
    _setIndexBase(_firstIndex);
    ASSERT(_data);
#ifndef NDEBUG
    ASSERT(p==data());
#endif
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, UpLo, Order, I, A>::_allocate(const ElementType &value)
{
    const IndexType numElements = 0.5*(_dim+1)*_dim;
    
    if (numElements==0) {
        return;
    }

    _raw_allocate();
    T *p = data();
    for (IndexType i=0; i<numElements; ++i) {
        _allocator.construct(p++, value);
    }
}

template <typename T, StorageUpLo UpLo, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, UpLo, Order, I, A>::_release()
{
    if (_data) {
        T *p = data();
        for (IndexType i=0; 2*i<(_dim+1)*_dim; ++i) {
            _allocator.destroy(p++);
        }
        _allocator.deallocate(data(), 0.5*((_dim+1)*_dim));
        _data = 0;
    }
    ASSERT(_data==0);
}

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_TCC