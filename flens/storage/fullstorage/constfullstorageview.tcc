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

namespace flens {

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(
                                        IndexType numRows, IndexType numCols,
                                        const ElementType *data,
                                        IndexType leadingDimension,
                                        IndexType firstRow,
                                        IndexType firstCol,
                                        const Allocator &allocator)
    : _data(data), _allocator(allocator), _numRows(numRows), _numCols(numCols),
      _leadingDimension(leadingDimension), _firstRow(0), _firstCol(0)
{
    changeIndexBase(firstRow, firstCol);
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(
                                                const ConstFullStorageView &rhs)
    : _data(rhs._data),
      _allocator(rhs._allocator),
      _numRows(rhs._numRows), _numCols(rhs._numCols),
      _leadingDimension(rhs._leadingDimension),
      _firstRow(rhs._firstRow), _firstCol(rhs._firstCol)
{
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
template <typename RHS>
ConstFullStorageView<T, Order, I, A>::ConstFullStorageView(const RHS &rhs)
    : _data(rhs.data()),
      _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _leadingDimension(rhs.leadingDimension()),
      _firstRow(0), _firstCol(0)
{
    changeIndexBase(rhs.firstRow(), rhs.firstCol());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
ConstFullStorageView<T, Order, I, A>::~ConstFullStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ElementType &
ConstFullStorageView<T, Order, I, A>::operator()(IndexType row,
                                                 IndexType col) const
{
    ASSERT(row>=_firstRow);
    ASSERT(row<_firstRow+_numRows);
    ASSERT(col>=_firstCol);
    ASSERT(col<_firstCol+_numCols);

    if (Order==cxxblas::ColMajor) {
        return _data[col*_leadingDimension+row];
    }
    return _data[row*_leadingDimension+col];
}

//-- methods -------------------------------------------------------------------

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::firstRow() const
{
    return _firstRow;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::firstCol() const
{
    return _firstCol;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::numRows() const
{
    return _numRows;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::numCols() const
{
    return _numCols;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::leadingDimension() const
{
    return _leadingDimension;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::strideRow() const
{
    return (Order==cxxblas::ColMajor) ? 1
                                      : leadingDimension();
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
typename ConstFullStorageView<T, Order, I, A>::IndexType
ConstFullStorageView<T, Order, I, A>::strideCol() const
{
    return (Order==cxxblas::ColMajor) ? leadingDimension()
                                      : 1;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ElementType *
ConstFullStorageView<T, Order, I, A>::data() const
{
#   ifndef NDEBUG
    if ((numRows()==0) || numCols()==0) {
        return 0;
    }
#   endif

    return &(this->operator()(_firstRow, _firstCol));
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::Allocator &
ConstFullStorageView<T, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
void
ConstFullStorageView<T, Order, I, A>::changeIndexBase(IndexType firstRow,
                                                      IndexType firstCol)
{
    if (Order==cxxblas::RowMajor) {
        _data = data() - (firstRow*leadingDimension() + firstCol);
    }
    if (Order==cxxblas::ColMajor) {
        _data = data() - (firstCol*leadingDimension() + firstRow);
    }
    _firstRow = firstRow;
    _firstCol = firstCol;
}

// view of rectangular part
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const ConstFullStorageView<T, Order, I, A>
ConstFullStorageView<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                           IndexType toRow, IndexType toCol,
                                           IndexType firstViewRow,
                                           IndexType firstViewCol) const
{
    const IndexType numRows = toRow-fromRow+1;
    const IndexType numCols = toCol-fromCol+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if ((numRows==0) || (numCols==0)) {
        fromRow = toRow = this->firstRow();
        fromCol = toCol = this->firstCol();
    }
#   endif

    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    return ConstView(numRows,                              // # rows
                     numCols,                              // # cols
                     &(this->operator()(fromRow, fromCol)),// data
                     leadingDimension(),                   // leading dimension
                     firstViewRow,                         // firstRow
                     firstViewCol,                         // firstCol
                     allocator());                         // allocator
}

// view of single row
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                              IndexType firstViewIndex) const
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numCols()==0) {
        row = this->firstRow();
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(numCols(),
                          &(this->operator()(row, _firstCol))-firstViewIndex,
                          strideCol(),
                          firstViewIndex,
                          allocator());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                              IndexType firstCol,
                                              IndexType lastCol,
                                              IndexType firstViewIndex) const
{
    const IndexType length = lastCol-firstCol+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        row =       firstRow();
        firstCol =  this->firstCol();
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ConstArrayView(length,
                          &(this->operator()(row, firstCol)),
                          strideCol(),
                          firstViewIndex,
                          allocator());
}

// view of single columns
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewCol(IndexType col,
                                              IndexType firstViewIndex) const
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numRows()==0) {
        col = this->firstCol();
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(numRows(),
                          &(this->operator()(_firstRow, col))-firstViewIndex,
                          strideRow(),
                          firstViewIndex,
                          allocator());
}

template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewCol(IndexType firstRow,
                                              IndexType lastRow,
                                              IndexType col,
                                              IndexType firstViewIndex) const
{
    const IndexType length = lastRow-firstRow+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        firstRow =  this->firstRow();
        col =       firstCol();
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ConstArrayView(length,
                          &(this->operator()(firstRow, col)),
                          strideRow(),
                          firstViewIndex,
                          allocator());
}

// view of d-th diagonal
template <typename T, cxxblas::StorageOrder Order, typename I, typename A>
const typename ConstFullStorageView<T, Order, I, A>::ConstArrayView
ConstFullStorageView<T, Order, I, A>::viewDiag(IndexType d,
                                               IndexType firstViewIndex) const
{
    IndexType col = firstCol() + ( (d>0) ? d : 0 );
    IndexType row = firstRow() + ( (d>0) ? 0 : -d );

    return ConstArrayView(std::min(numRows(),numCols()) - std::abs(d),
                          &(this->operator()(row,col)),
                          leadingDimension()+1,
                          firstViewIndex,
                          allocator());
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_CONSTFULLSTORAGEVIEW_TCC
