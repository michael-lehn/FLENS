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

#ifndef FLENS_STORAGE_FULLSTORAGE_FULLSTORAGEVIEW_TCC
#define FLENS_STORAGE_FULLSTORAGE_FULLSTORAGEVIEW_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/fullstorage/fullstorageview.h>
#include <flens/storage/fullstorage/trapezoidalfill.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, StorageOrder Order, typename I, typename A>
FullStorageView<T, Order, I, A>::FullStorageView(IndexType numRows,
                                                 IndexType numCols,
                                                 ElementType *data,
                                                 IndexType leadingDimension,
                                                 IndexType firstRow,
                                                 IndexType firstCol,
                                                 const Allocator &allocator)
    : data_(data),
      allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      leadingDimension_(leadingDimension),
      firstRow_(0), firstCol_(0)
{
    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);

    changeIndexBase(firstRow, firstCol);
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename ARRAY>
FullStorageView<T, Order, I, A>::FullStorageView(IndexType numRows,
                                                 IndexType numCols,
                                                 ARRAY &array,
                                                 IndexType leadingDimension,
                                                 IndexType firstRow,
                                                 IndexType firstCol,
                                                 const Allocator &allocator)
    : data_(array.data()),
      allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      leadingDimension_(leadingDimension),
      firstRow_(0), firstCol_(0)
{
    ASSERT(numRows*numCols<=array.length());
    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);

    changeIndexBase(firstRow, firstCol);
}

template <typename T, StorageOrder Order, typename I, typename A>
FullStorageView<T, Order, I, A>::FullStorageView(const FullStorageView &rhs)
    : data_(rhs.data_),
      allocator_(rhs.allocator_),
      numRows_(rhs.numRows_), numCols_(rhs.numCols_),
      leadingDimension_(rhs.leadingDimension_),
      firstRow_(rhs.firstRow_), firstCol_(rhs.firstCol_)
{
    ASSERT(order==rhs.order);
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
FullStorageView<T, Order, I, A>::FullStorageView(RHS &rhs)
    : data_(rhs.data()),
      allocator_(rhs.allocator()),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      leadingDimension_(rhs.leadingDimension()),
      firstRow_(0), firstCol_(0)
{
    ASSERT(order==rhs.order);
    changeIndexBase(rhs.firstRow(), rhs.firstCol());
}

template <typename T, StorageOrder Order, typename I, typename A>
FullStorageView<T, Order, I, A>::~FullStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ElementType &
FullStorageView<T, Order, I, A>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (numRows()>0 && numCols()>0) {
        ASSERT(row>=firstRow_);
        ASSERT(row<firstRow_+numRows_);
        ASSERT(col>=firstCol_);
        ASSERT(col<firstCol_+numCols_);
        ASSERT(data_);
    } else {
        ASSERT(row==firstRow_);
        ASSERT(col==firstCol_);
    }
#   endif

    if (Order==ColMajor) {
        return data_[col*leadingDimension_+row];
    }
    return data_[row*leadingDimension_+col];
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ElementType &
FullStorageView<T, Order, I, A>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (numRows()>0 && numCols()>0) {
        ASSERT(row>=firstRow_);
        ASSERT(row<firstRow_+numRows_);
        ASSERT(col>=firstCol_);
        ASSERT(col<firstCol_+numCols_);
        ASSERT(data_);
    } else {
        ASSERT(row==firstRow_);
        ASSERT(col==firstCol_);
    }
#   endif

    if (Order==ColMajor) {
        return data_[col*leadingDimension_+row];
    }
    return data_[row*leadingDimension_+col];
}

//-- Methods -------------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::firstRow() const
{
    return firstRow_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::lastRow() const
{
    return firstRow_+numRows_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::firstCol() const
{
    return firstCol_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::lastCol() const
{
    return firstCol_+numCols_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::numRows() const
{
    return numRows_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::numCols() const
{
    return numCols_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::leadingDimension() const
{
    return leadingDimension_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension();
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::IndexType
FullStorageView<T, Order, I, A>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()
                             : 1;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ElementType *
FullStorageView<T, Order, I, A>::data() const
{
    return &(operator()(firstRow_, firstCol_));
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ElementType *
FullStorageView<T, Order, I, A>::data()
{
    return &(operator()(firstRow_, firstCol_));
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::Allocator &
FullStorageView<T, Order, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
FullStorageView<T, Order, I, A>::resize(IndexType DEBUG_VAR(numRows_),
                                        IndexType DEBUG_VAR(numCols_),
                                        IndexType firstRow,
                                        IndexType firstCol,
                                        const ElementType &)
{
    ASSERT(numRows_==numRows());
    ASSERT(numCols_==numCols());

    changeIndexBase(firstRow, firstCol);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename FS>
bool
FullStorageView<T, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.numRows(), rhs.numCols(),
                  rhs.firstRow(), rhs.firstCol(),
                  value);
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
FullStorageView<T, Order, I, A>::fill(const ElementType &value)
{
    if (Order==RowMajor) {
        ElementType *p = data();
        for (IndexType i=0; i<numRows(); ++i, p+=leadingDimension()) {
            std::fill_n(p, numCols(), value);
        }
        return true;
    }
    if (Order==ColMajor) {
        ElementType *p = data();
        for (IndexType j=0; j<numCols(); ++j, p+=leadingDimension()) {
            std::fill_n(p, numRows(), value);
        }
        return true;
    }
    ASSERT(0);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
FullStorageView<T, Order, I, A>::fill(StorageUpLo  upLo,
                                      const ElementType &value)
{
    trapezoidalFill(order, upLo, value,
                    numRows(), numCols(),
                    data(), leadingDimension());
    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
FullStorageView<T, Order, I, A>::changeIndexBase(IndexType firstRow,
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
        if (Order==RowMajor) {
            data_ = data() - (firstRow*leadingDimension() + firstCol);
        }
        if (Order==ColMajor) {
            data_ = data() - (firstCol*leadingDimension() + firstRow);
        }
    }
    firstRow_ = firstRow;
    firstCol_ = firstCol;
}

// view of fullstorage scheme as an array
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::arrayView(IndexType firstViewIndex) const
{
#   ifndef NDEBUG
    if (order==ColMajor) {
        ASSERT(numRows()==leadingDimension());
    } else {
        ASSERT(numCols()==leadingDimension());
    }
#   endif

    return ConstArrayView(numCols()*numRows(),
                          data(),
                          IndexType(1),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::arrayView(IndexType firstViewIndex)
{
#   ifndef NDEBUG
    if (order==ColMajor) {
        ASSERT(numRows()==leadingDimension());
    } else {
        ASSERT(numCols()==leadingDimension());
    }
#   endif

    return ArrayView(numCols()*numRows(),
                     data(),
                     IndexType(1),
                     firstViewIndex,
                     allocator());
}


// view of rectangular part
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstView
FullStorageView<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                      IndexType toRow, IndexType toCol,
                                      IndexType strideRow, IndexType strideCol,
                                      IndexType firstViewRow,
                                      IndexType firstViewCol) const
{
    const IndexType numRows = (toRow-fromRow)/strideRow+1;
    const IndexType numCols = (toCol-fromCol)/strideCol+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if ((numRows==0) || (numCols==0)) {
        return ConstView(numRows, numCols, 0, leadingDimension(),
                         firstViewRow, firstViewCol, allocator());
    }

#   endif

    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    ASSERT(order==ColMajor || strideCol==IndexType(1) );
    ASSERT(order==RowMajor || strideRow==IndexType(1) );

    return ConstView(numRows,                               // # rows
                     numCols,                               // # cols
                     &(operator()(fromRow, fromCol)),       // data
                     leadingDimension()*strideRow*strideCol,// leading dimension
                     firstViewRow,                          // firstRow
                     firstViewCol,                          // firstCol
                     allocator());                          // allocator
}

// view of rectangular part
template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::View
FullStorageView<T, Order, I, A>::view(IndexType fromRow, IndexType fromCol,
                                      IndexType toRow, IndexType toCol,
                                      IndexType strideRow, IndexType strideCol,
                                      IndexType firstViewRow,
                                      IndexType firstViewCol)
{
    const IndexType numRows = (toRow-fromRow)/strideRow+1;
    const IndexType numCols = (toCol-fromCol)/strideCol+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if ((numRows==0) || (numCols==0)) {
        return      View(numRows, numCols, 0, leadingDimension(),
                         firstViewRow, firstViewCol, allocator());
    }

#   endif

    ASSERT(fromRow>=firstRow());
    ASSERT(fromRow<=toRow);
    ASSERT(toRow<=lastRow());

    ASSERT(fromCol>=firstCol());
    ASSERT(fromCol<=toCol);
    ASSERT(toCol<=lastCol());

    ASSERT(order==ColMajor || strideCol==IndexType(1) );
    ASSERT(order==RowMajor || strideRow==IndexType(1) );

    return      View(numRows,                               // # rows
                     numCols,                               // # cols
                     &(operator()(fromRow, fromCol)),       // data
                     leadingDimension()*strideRow*strideCol,// leading dimension
                     firstViewRow,                          // firstRow
                     firstViewCol,                          // firstCol
                     allocator());                          // allocator
}

// view of single row
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewRow(IndexType row,
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
                          &(operator()(row, firstCol_)),
                          strideCol(),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                         IndexType firstViewIndex)
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numCols()==0) {
        return ArrayView(numCols(), 0, strideCol(),
                         firstViewIndex, allocator());
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ArrayView(numCols(),
                     &(operator()(row, firstCol_)),
                     strideCol(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                         IndexType firstCol, IndexType lastCol,
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
                          &(operator()(row, firstCol)),
                          strideCol()*stride,
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewRow(IndexType row,
                                         IndexType firstCol, IndexType lastCol,
                                         IndexType stride,
                                         IndexType firstViewIndex)
{
    const IndexType length = (lastCol-firstCol)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ArrayView(length, 0, strideCol()*stride,
                         firstViewIndex, allocator());
    }
#   endif

    ASSERT(row>=firstRow());
    ASSERT(row<=lastRow());

    return ArrayView(length,
                     &(operator()(row, firstCol)),
                     strideCol()*stride,
                     firstViewIndex,
                     allocator());
}

// view of single column
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewCol(IndexType col,
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
                          &(operator()(firstRow_, col)),
                          strideRow(),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewCol(IndexType col,
                                         IndexType firstViewIndex)
{
#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (numRows()==0) {
        return ArrayView(numRows(), 0, strideRow(),
                         firstViewIndex, allocator());
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ArrayView(numRows(),
                     &(operator()(firstRow_, col)),
                     strideRow(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
                                         IndexType stride, IndexType col,
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
                          &(operator()(firstRow, col)),
                          strideRow()*stride,
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
                                         IndexType stride, IndexType col,
                                         IndexType firstViewIndex)
{
    const IndexType length = (lastRow-firstRow)/stride+1;

#   ifndef NDEBUG
    // prevent an out-of-bound assertion in case a view is empty anyway
    if (length==0) {
        return ArrayView(length, 0, strideRow()*stride,
                         firstViewIndex, allocator());
    }
#   endif

    ASSERT(col>=firstCol());
    ASSERT(col<=lastCol());

    return ArrayView(length,
                     &(operator()(firstRow, col)),
                     strideRow()*stride,
                     firstViewIndex,
                     allocator());
}

// view of d-th diagonal
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewDiag(IndexType d,
                                          IndexType firstViewIndex) const
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ConstArrayView(std::min(numRows()-row_, numCols()-col_),
                          &(operator()(row,col)),
                          leadingDimension()+1,
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewDiag(IndexType d,
                                          IndexType firstViewIndex)
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ArrayView(std::min(numRows()-row_,numCols()-col_),
                     &(operator()(row,col)),
                     leadingDimension()+1,
                     firstViewIndex,
                     allocator());
}

// view of d-th anti-diagonal
template <typename T, StorageOrder Order, typename I, typename A>
const typename FullStorageView<T, Order, I, A>::ConstArrayView
FullStorageView<T, Order, I, A>::viewAntiDiag(IndexType d,
                                             IndexType firstViewIndex) const
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ConstArrayView(std::min(numRows()-row_, numCols()-col_),
                          &(operator()(row,lastCol()-col+1)),
                          -leadingDimension()+1,
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename FullStorageView<T, Order, I, A>::ArrayView
FullStorageView<T, Order, I, A>::viewAntiDiag(IndexType d,
                                              IndexType firstViewIndex)
{
    IndexType row_ = (d>0) ? 0 : -d;
    IndexType col_ = (d>0) ? d :  0;

    IndexType row = firstRow() + row_;
    IndexType col = firstCol() + col_;

    return ArrayView(std::min(numRows()-row_,numCols()-col_),
                     &(operator()(row,lastCol()-col+1)),
                     -leadingDimension()+1,
                     firstViewIndex,
                     allocator());
}


//-- FullStorageView specific functions ----------------------------------------

//
//  fillRandom
//

template <typename T, StorageOrder Order, typename I, typename Allocator>
bool
fillRandom(FullStorageView<T, Order, I, Allocator> &A)
{
    typedef FullStorageView<T,Order,I,Allocator>    FullStorageView;
    typedef typename FullStorageView::ElementType   ElementType;
    typedef typename FullStorageView::IndexType     IndexType;

    if (Order==RowMajor) {
        ElementType *p = A.data();
        for (IndexType i=0; i<A.numRows(); ++i, p+=A.leadingDimension()) {
            for (IndexType j=0;j<A.numCols();++j) {
                p[j] = randomValue<T>();
            }
        }
        return true;
    }
    if (Order==ColMajor) {
        ElementType *p = A.data();
        for (IndexType j=0; j<A.numCols(); ++j, p+=A.leadingDimension()) {
            for (IndexType i=0;i<A.numRows();++i) {
                p[i] = randomValue<T>();
            }
        }
        return true;
    }
    ASSERT(0);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename Allocator>
bool
fillRandom(StorageUpLo upLo, FullStorageView<T, Order, I, Allocator> &A)
{
    trapezoidalFillRandom(Order, upLo,
                          A.numRows(), A.numCols(),
                          A.data(), A.leadingDimension());
    return true;
}

} // namespace flens

#endif // FLENS_STORAGE_FULLSTORAGE_FULLSTORAGEVIEW_TCC
