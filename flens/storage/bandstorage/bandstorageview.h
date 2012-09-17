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

#ifndef FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_H
#define FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_H 1

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

template <typename T, StorageOrder Order, typename I, typename A>
    class BandStorage;

template <typename T, StorageOrder Order, typename I, typename A>
    class ConstBandStorageView;

template <typename T,
          StorageOrder Order = ColMajor,
          typename I = IndexOptions<>,
          typename A = std::allocator<T> >
class BandStorageView
{
    public:
        typedef T                                     ElementType;
        typedef typename I::IndexType                 IndexType;
        typedef A                                     Allocator;

        static const StorageOrder                     order = Order;
        static const IndexType                        defaultIndexBase
                                                          = I::defaultIndexBase;

        typedef ConstBandStorageView<T, Order, I, A>  ConstView;
        typedef BandStorageView                       View;
        typedef BandStorage<T, Order, I, A>           NoView;

        typedef flens::ConstArrayView<T, I, A>        ConstArrayView;
        typedef flens::ArrayView<T, I, A>             ArrayView;
        typedef flens::Array<T, I, A>                 Array;

        BandStorageView(IndexType numRows, IndexType numCols,
                        IndexType numSubDiags, IndexType numSuperDiags,
                        ElementType *data,
                        IndexType leadingDimension,
                        IndexType firstIndex = I::defaultIndexBase,
                        const Allocator &allocator = Allocator());


        BandStorageView(const BandStorageView &rhs);

        template <typename RHS>
            BandStorageView(RHS &rhs);

        ~BandStorageView();

        //-- operators ---------------------------------------------------------

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        //-- methods -----------------------------------------------------------

        IndexType
        firstRow() const;

        IndexType
        lastRow() const;

        IndexType
        firstCol() const;

        IndexType
        lastCol() const;

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        IndexType
        numRows() const;

        IndexType
        numCols() const;

        IndexType
        dim() const;

        IndexType
        numSubDiags() const;

        IndexType
        numSuperDiags() const;

        IndexType
        leadingDimension() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        const Allocator &
        allocator() const;

        bool
        resize(IndexType numRows, IndexType numCols,
               IndexType numSubDiags, IndexType numSuperDiags,
               IndexType firstIndex = I::defaultIndexBase,
               const ElementType &value = ElementType());

        template <typename FS>
            bool
            resize(const FS &rhs, const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        bool
        fillRandom();

        void
        changeIndexBase(IndexType firstIndex);

        // view as an array
        const ConstArrayView
        arrayView(IndexType firstViewIndex = I::defaultIndexBase) const ;

        ArrayView
        arrayView(IndexType firstViewIndex = I::defaultIndexBase);

        // view of a diagonal
        const ConstArrayView
        viewDiag(IndexType diag,
                 IndexType firstViewIndex = I::defaultIndexBase) const;

        ArrayView
        viewDiag(IndexType diag,
                 IndexType firstViewIndex = I::defaultIndexBase);

        // view of some diagonals
        const ConstView
        viewDiags(IndexType fromDiag, IndexType toDiag) const;

        View
        viewDiags(IndexType fromDiag, IndexType toDiag);

    private:
        ElementType  *_data;
        Allocator    _allocator;
        IndexType    _numRows, _numCols;
        IndexType    _numSubDiags, _numSuperDiags;
        IndexType    _firstIndex;
        IndexType    _leadingDimension;
};

} // namespace flens

#endif // FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_H
