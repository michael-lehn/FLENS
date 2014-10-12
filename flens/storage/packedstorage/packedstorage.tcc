/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
view *
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

#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/packedstorage/packedstorage.h>
#include <flens/typedefs.h>

namespace flens {

//-- Constructors --------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorage<T, Order, I, A>::PackedStorage(IndexType dim,
                                             IndexType indexBase,
                                             const ElementType &value,
                                             const Allocator &allocator)
    : data_(0), allocator_(allocator),
      dim_(dim),
      indexBase_(indexBase)
{
    ASSERT(dim_>=0);

    allocate_(value);
}

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorage<T, Order, I, A>::PackedStorage(const PackedStorage &rhs)
    : data_(0), allocator_(rhs.allocator()),
      dim_(rhs.dim()),
      indexBase_(rhs.indexBase())
{
    allocate_(ElementType());

    cxxblas::copy(rhs.dim(), rhs.data(), 1, data(), 1);
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
PackedStorage<T, Order, I, A>::PackedStorage(const RHS &rhs)
    : data_(0), allocator_(rhs.allocator()),
      dim_(rhs.dim()),
      indexBase_(rhs.indexBase())
{
    allocate_(ElementType());

    cxxblas::copy(rhs.dim(), rhs.data(), 1, data(), 1);

}

template <typename T, StorageOrder Order, typename I, typename A>
PackedStorage<T, Order, I, A>::~PackedStorage()
{
    release_();
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, Order, I, A>::ElementType &
PackedStorage<T, Order, I, A>::operator()(StorageUpLo  upLo,
                                          IndexType    row,
                                          IndexType    col) const
{
    ASSERT(row>=indexBase());
    ASSERT(col>=indexBase());
    ASSERT(row<=indexBase()+dim()-1);
    ASSERT(col<=indexBase()+dim()-1);

    if (upLo==Lower) {
        ASSERT(row>=col);
    } else {
        ASSERT(col>=row);
    }

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

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, Order, I, A>::ElementType &
PackedStorage<T, Order, I, A>::operator()(StorageUpLo  upLo,
                                          IndexType    row,
                                          IndexType    col)
{
    ASSERT(row>=indexBase());
    ASSERT(col>=indexBase());
    ASSERT(row<=indexBase()+dim()-1);
    ASSERT(col<=indexBase()+dim()-1);

    if (upLo==Lower) {
        ASSERT(row>=col);
    } else {
        ASSERT(col>=row);
    }

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

//-- Methods -------------------------------------------------------------------
template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, Order, I, A>::IndexType
PackedStorage<T, Order, I, A>::indexBase() const
{
    return indexBase_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, Order, I, A>::IndexType
PackedStorage<T, Order, I, A>::numNonZeros() const
{
    return (dim_+1)*dim_/IndexType(2);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, Order, I, A>::IndexType
PackedStorage<T, Order, I, A>::dim() const
{
    return dim_;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, Order, I, A>::ElementType *
PackedStorage<T, Order, I, A>::data() const
{
    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename PackedStorage<T, Order, I, A>::ElementType *
PackedStorage<T, Order, I, A>::data()
{
    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename PackedStorage<T, Order, I, A>::Allocator &
PackedStorage<T, Order, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
PackedStorage<T, Order, I, A>::resize(IndexType dim,
                                            IndexType indexBase,
                                            const ElementType &value)
{
    if (dim_!=dim) {
        release_();
        dim_ = dim;
        indexBase_ = indexBase;
        allocate_(value);
        return true;
    }
    changeIndexBase(indexBase);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
PackedStorage<T, Order, I, A>::fill(const ElementType &value)
{
    ASSERT(data_);
    std::fill_n(data(), 0.5*(dim()+1)*dim(), value);
    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, Order, I, A>::changeIndexBase(IndexType indexBase)
{
    indexBase_ = indexBase;
}

//-- Private Methods -----------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, Order, I, A>::raw_allocate_()
{
    ASSERT(!data_);
    ASSERT(dim_>0);

    data_ = allocator_.allocate(numNonZeros());
    ASSERT(data_);

#ifndef NDEBUG
    ElementType *p = data_;
#endif

    changeIndexBase(indexBase_);

#ifndef NDEBUG
    ASSERT(p==data());
#endif
}

template <typename T, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, Order, I, A>::allocate_(const ElementType &value)
{
    const IndexType numElements = numNonZeros();

    if (numElements==0) {
        return;
    }

    raw_allocate_();
    T *p = data();
    for (IndexType i=0; i<numElements; ++i, ++p) {
        allocator_.construct(p, value);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
void
PackedStorage<T, Order, I, A>::release_()
{
    if (data_) {
        const IndexType numElements = numNonZeros();

        T *p = data();
        for (IndexType i=0; i<numElements; ++i, ++p) {
            allocator_.destroy(p);
        }
        allocator_.deallocate(data(), numElements);
        data_ = 0;
    }
    ASSERT(data_==0);
}

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_TCC
