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

#ifndef FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_H
#define FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_H 1

#include <cxxblas/typedefs.h>
#include <flens/typedefs.h>
#include <flens/storage/indexoptions.h>

namespace flens {

template <typename T, typename I, typename A>
    class Array;

template <typename T, typename I, typename A>
    class ArrayView;

template <typename T, typename I, typename A>
    class ConstArrayView;

template <typename T, StorageUpLo UpLo, StorageOrder Order,
          typename I, typename A>
    class PackedStorage;

template <typename T, StorageUpLo UpLo, StorageOrder Order,
          typename I, typename A>
    class PackedStorageView;

template <typename T,
          StorageUpLo UpLo = Upper,
          StorageOrder Order = ColMajor,
          typename I = IndexOptions<>,
          typename A = std::allocator<T> >
class ConstPackedStorageView
{
    public:
        typedef T                        ElementType;
        typedef typename I::IndexType    IndexType;
        typedef A                        Allocator;

        static const StorageOrder        order = Order;
        static const StorageUpLo         upLo  = UpLo;
        static const IndexType           defaultIndexBase = I::defaultIndexBase;

        typedef ConstPackedStorageView                   ConstView;
        typedef PackedStorageView<T, UpLo, Order, I, A>  View;
        typedef PackedStorage<T, UpLo, Order, I, A>      NoView;

        typedef flens::ConstArrayView<T, I, A>           ConstArrayView;
        typedef flens::ArrayView<T, I, A>                ArrayView;
        typedef flens::Array<T, I, A>                    Array;

        ConstPackedStorageView(IndexType dim,
                               const ElementType *data,
                               IndexType indexBase = I::defaultIndexBase,
                               const Allocator &allocator = Allocator());

        template <typename ARRAY>
            ConstPackedStorageView(IndexType dim,
                                   ARRAY &array,
                                   IndexType indexBase = I::defaultIndexBase,
                                   const Allocator &allocator = Allocator());

        ConstPackedStorageView(const ConstPackedStorageView &rhs);

        template <typename RHS>
            ConstPackedStorageView(const RHS &rhs);

        ~ConstPackedStorageView();

        //-- operators ---------------------------------------------------------

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        //-- methods -----------------------------------------------------------

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        IndexType
        dim() const;

        const ElementType *
        data() const;

        const Allocator &
        allocator() const;

        void
        changeIndexBase(IndexType indexBase);


    private:
        const ElementType    *_data;
        Allocator            _allocator;
        IndexType            _dim;
        IndexType            _firstIndex;
};

} // namespace flens

#endif // FLENS_STORAGE_PACKEDSTORAGE_CONSTPACKEDSTORAGEVIEW_H
