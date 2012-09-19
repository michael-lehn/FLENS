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

#ifndef FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_H
#define FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_H

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/indexoptions.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, typename I, typename A>
    class Array;

template <typename T, typename I, typename A>
    class ArrayView;

template <typename T, typename I, typename A>
    class ConstArrayView;

template <typename T, StorageOrder Order, typename I, typename A>
    class PackedStorageView;

template <typename T, StorageOrder Order, typename I, typename A>
    class ConstPackedStorageView;

template <typename T,
          StorageOrder Order = ColMajor,
          typename I = IndexOptions<>,
          typename A = std::allocator<T> >
class PackedStorage
{
    public:
        typedef T                        ElementType;
        typedef typename I::IndexType    IndexType;
        typedef A                        Allocator;

        static const StorageOrder        order = Order;
        static const IndexType           defaultIndexBase = I::defaultIndexBase;

        typedef ConstPackedStorageView<T, Order, I, A>  ConstView;
        typedef PackedStorageView<T, Order, I, A>       View;
        typedef PackedStorage                           NoView;

        typedef flens::ConstArrayView<T, I, A>          ConstArrayView;
        typedef flens::ArrayView<T, I, A>               ArrayView;
        typedef flens::Array<T, I, A>                   Array;

        PackedStorage(IndexType dim,
                      IndexType indexBase = I::defaultIndexBase,
                      const ElementType &value = ElementType(),
                      const Allocator &allocator = Allocator());

        PackedStorage(const PackedStorage &rhs);

        template <typename RHS>
            PackedStorage(const RHS &rhs);

        ~PackedStorage();

        //-- operators ---------------------------------------------------------

        const ElementType &
        operator()(StorageUpLo  upLo, IndexType row, IndexType col) const;

        ElementType &
        operator()(StorageUpLo  upLo, IndexType row, IndexType col);

        //-- methods -----------------------------------------------------------


        IndexType
        indexBase() const;

        IndexType
        numNonZeros() const;

        IndexType
        dim() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        const Allocator &
        allocator() const;

        bool
        resize(IndexType dim,
               IndexType indexBase = I::defaultIndexBase,
               const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        void
        changeIndexBase(IndexType indexBase);


    private:

        void
        _raw_allocate();

        void
        _allocate(const ElementType &value = ElementType());

        void
        _release();

        ElementType  *_data;
        Allocator    _allocator;
        IndexType    _dim;
        IndexType    _indexBase;
};

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_PACKEDSTORAGE_H
