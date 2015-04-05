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

#ifndef FLENS_STORAGE_FULLSTORAGE_CONSTFULLSTORAGEVIEW_TCC
#define FLENS_STORAGE_FULLSTORAGE_CONSTFULLSTORAGEVIEW_TCC 1

#include <flens/storage/fullstorage/constfullstorageview.h>
#include <flens/auxiliary/macros.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, StorageOrder NVOrder, typename I, typename A>
ConstFullStorageView<T, NVOrder, I, A>::ConstFullStorageView(
                                                     IndexType numRows,
                                                     IndexType numCols,
                                                     IndexType leadingDimension,
                                                     const ElementType *data,
                                                     IndexType firstRow,
                                                     IndexType firstCol,
                                                     const Allocator &allocator)
    : data_(data), allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      strideRow_((NVOrder==ColMajor) ? IndexType(1) : leadingDimension),
      strideCol_((NVOrder==ColMajor) ? leadingDimension : IndexType(1)),
      firstRow_(0), firstCol_(0)
{
#   ifndef NDEBUG
    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);
    if (NVOrder==ColMajor) {
        ASSERT(leadingDimension>=numRows);
    }
    if (NVOrder==RowMajor) {
        ASSERT(leadingDimension>=numCols);
    }
#   endif

    changeIndexBase(firstRow, firstCol);
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
ConstFullStorageView<T, NVOrder, I, A>::ConstFullStorageView(
                                                     IndexType numRows,
                                                     IndexType numCols,
                                                     const ElementType *data,
                                                     IndexType strideRow,
                                                     IndexType strideCol,
                                                     IndexType firstRow,
                                                     IndexType firstCol,
                                                     const Allocator &allocator)
    : data_(data), allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      strideRow_(strideRow), strideCol_(strideCol),
      firstRow_(0), firstCol_(0)
{
    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);
    changeIndexBase(firstRow, firstCol);
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
ConstFullStorageView<T, NVOrder, I, A>::ConstFullStorageView(
                                                const ConstFullStorageView &rhs)
    : data_(rhs.data_),
      allocator_(rhs.allocator_),
      numRows_(rhs.numRows_), numCols_(rhs.numCols_),
      strideRow_(rhs.strideRow()), strideCol_(rhs.strideCol()),
      firstRow_(rhs.firstRow_), firstCol_(rhs.firstCol_)
{
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
template <typename RHS>
ConstFullStorageView<T, NVOrder, I, A>::ConstFullStorageView(const RHS &rhs)
    : data_(rhs.data()),
      allocator_(rhs.allocator()),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      strideRow_(rhs.strideRow()), strideCol_(rhs.strideCol()),
      firstRow_(0), firstCol_(0)
{
    changeIndexBase(rhs.firstRow(), rhs.firstCol());
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
ConstFullStorageView<T, NVOrder, I, A>::~ConstFullStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ElementType &
ConstFullStorageView<T, NVOrder, I, A>::operator()(IndexType row,
                                                   IndexType col) const
{
#   ifndef NDEBUG
    if (numRows()>0 && numCols()>0) {
        ASSERT(row>=firstRow_);
        ASSERT(row<firstRow_+numRows_);
        ASSERT(col>=firstCol_);
        ASSERT(col<firstCol_+numCols_);
    } else {
        ASSERT(row==firstRow_);
        ASSERT(col==firstCol_);
        ASSERT(data_);
    }
#   endif

    if (strideRow_==IndexType(1)) {
        return data_[row+col*strideCol_];
    }
    if (strideCol_==IndexType(1)) {
        return data_[row*strideRow_+col];
    }
    return data_[row*strideRow_+col*strideCol_];
}

//-- methods -------------------------------------------------------------------

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::firstRow() const
{
    return firstRow_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::firstCol() const
{
    return firstCol_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::lastRow() const
{
    return firstRow_+numRows_-1;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::lastCol() const
{
    return firstCol_+numCols_-1;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::numRows() const
{
    return numRows_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::numCols() const
{
    return numCols_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
StorageOrder
ConstFullStorageView<T, NVOrder, I, A>::order() const
{
    if (strideRow_==IndexType(1)) {
        return ColMajor;
    }
    if (strideCol_==IndexType(1)) {
        return RowMajor;
    }
    return Grid;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::leadingDimension() const
{
    if (strideRow_==IndexType(1)) {
        return strideCol_;
    }
    if (strideCol_==IndexType(1)) {
        return strideRow_;
    }
    ASSERT(0); // Grid views do not have a leading dimension
    return 0;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::strideRow() const
{
    return strideRow_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
typename ConstFullStorageView<T, NVOrder, I, A>::IndexType
ConstFullStorageView<T, NVOrder, I, A>::strideCol() const
{
    return strideCol_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ElementType *
ConstFullStorageView<T, NVOrder, I, A>::data() const
{
    return &(this->operator()(firstRow_, firstCol_));
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::Allocator &
ConstFullStorageView<T, NVOrder, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
void
ConstFullStorageView<T, NVOrder, I, A>::changeIndexBase(IndexType firstRow,
                                                        IndexType firstCol)
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if ((numRows()==0) || (numCols()==0)) {
        firstRow_ = firstRow;
        firstCol_ = firstCol;
        return;
    }
#   endif

    if (data_) {
        data_ = data() - (firstRow*strideRow() + firstCol*strideCol());
    }
    firstRow_ = firstRow;
    firstCol_ = firstCol;
}

/// view of fullstorage scheme as an array
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::arrayView(IndexType firstViewIndex
                                                  ) const
{
#   ifndef NDEBUG
    if (numRows()*strideRow()!=strideCol()
     && numCols()*strideCol()!=strideRow())
    {
        ASSERT(0);
    }
#   endif

    return ConstArrayView(numCols()*numRows(),
                          data(),
                          std::min(strideRow(), strideCol()),
                          firstViewIndex,
                          allocator());
}

// view of rectangular part
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstView
ConstFullStorageView<T, NVOrder, I, A>::view(IndexType fromRow,
                                             IndexType fromCol,
                                             IndexType toRow,
                                             IndexType toCol,
                                             IndexType strideRow,
                                             IndexType strideCol,
                                             IndexType firstViewRow,
                                             IndexType firstViewCol) const
{
    const IndexType numRows = (toRow-fromRow)/strideRow+1;
    const IndexType numCols = (toCol-fromCol)/strideCol+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if ((numRows==0) || (numCols==0)) {
        return ConstView(numRows, numCols, 0,
                         this->strideRow(), this->strideCol(),
                         firstViewRow, firstViewCol, allocator());
    }

#   endif

    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    return ConstView(numRows,                               // # rows
                     numCols,                               // # cols
                     &(operator()(fromRow, fromCol)),       // data
                     this->strideRow()*strideRow,           // strideRow
                     this->strideCol()*strideCol,           // strideCol
                     firstViewRow,                          // firstRow
                     firstViewCol,                          // firstCol
                     allocator());                          // allocator
}

// view of single row
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewRow(IndexType row,
                                                IndexType firstViewIndex) const
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numCols()==0) {
        return ConstArrayView(numCols(), 0, strideCol(),
                              firstViewIndex, allocator());
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(numCols(),
                          &(this->operator()(row, firstCol_)),
                          strideCol(),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewRow(IndexType row,
                                                IndexType firstCol,
                                                IndexType lastCol,
                                                IndexType stride,
                                                IndexType firstViewIndex) const
{
    const IndexType length = (lastCol-firstCol)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ConstArrayView(length, 0, strideCol()*stride,
                              firstViewIndex, allocator());
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(length,
                          &(this->operator()(row, firstCol)),
                          strideCol()*stride,
                          firstViewIndex,
                          allocator());
}

// view of single columns
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewCol(IndexType col,
                                                IndexType firstViewIndex) const
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numRows()==0) {
        return ConstArrayView(numRows(), 0, strideRow(),
                              firstViewIndex, allocator());
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(numRows(),
                          &(this->operator()(firstRow_, col)),
                          strideRow(),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewCol(IndexType firstRow,
                                                IndexType lastRow,
                                                IndexType stride,
                                                IndexType col,
                                                IndexType firstViewIndex) const
{
    const IndexType length = (lastRow-firstRow)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ConstArrayView(length, 0, strideRow()*stride,
                              firstViewIndex, allocator());
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(length,
                          &(this->operator()(firstRow, col)),
                          strideRow()*stride,
                          firstViewIndex,
                          allocator());
}

// view of d-th diagonal
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewDiag(IndexType d,
                                                 IndexType firstViewIndex) const
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ConstArrayView(std::min(numRows()-row_, numCols()-col_),
                          &(this->operator()(row,col)),
                          strideRow()+strideCol(),
                          firstViewIndex,
                          allocator());
}

// view of d-th anti-diagonal
template <typename T, StorageOrder NVOrder, typename I, typename A>
const typename ConstFullStorageView<T, NVOrder, I, A>::ConstArrayView
ConstFullStorageView<T, NVOrder, I, A>::viewAntiDiag(IndexType d,
                                                     IndexType firstViewIndex
                                                     ) const
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ConstArrayView(std::min(numRows()-row_, numCols()-col_),
                          &(this->operator()(row,lastCol()-col+1)),
                          -strideCol()+strideRow(),
                          firstViewIndex,
                          allocator());
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_CONSTFULLSTORAGEVIEW_TCC
