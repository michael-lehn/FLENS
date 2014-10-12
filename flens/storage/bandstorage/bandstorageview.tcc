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

#ifndef FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_TCC
#define FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/bandstorage/bandstorageview.h>
#include <flens/typedefs.h>

namespace flens {

template <typename T, StorageOrder Order, typename I, typename A>
BandStorageView<T, Order, I, A>::BandStorageView(IndexType numRows,
                                                 IndexType numCols,
                                                 IndexType numSubDiags,
                                                 IndexType numSuperDiags,
                                                 ElementType *data,
                                                 IndexType leadingDimension,
                                                 IndexType firstIndex,
                                                 const Allocator &allocator)
    : data_(data),
      allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      numSubDiags_(numSubDiags), numSuperDiags_(numSuperDiags),
      firstIndex_(firstIndex),
      leadingDimension_(leadingDimension)
{
    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);
    ASSERT(numSubDiags_>=0);
    ASSERT(numSuperDiags_>=0);

}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorageView<T, Order, I, A>::BandStorageView(const BandStorageView &rhs)
    : data_(rhs.data_),
      allocator_(rhs.allocator_),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      numSubDiags_(rhs.numSubDiags_), numSuperDiags_(rhs.numSuperDiags_),
      firstIndex_(rhs.firstIndex()),
      leadingDimension_(rhs.leadingDimension_)
{

}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
BandStorageView<T, Order, I, A>::BandStorageView(RHS &rhs)
    : data_(rhs.data()),
      allocator_(rhs.allocator()),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      numSubDiags_(rhs.numSubDiags()), numSuperDiags_(rhs.numSuperDiags()),
      firstIndex_(rhs.firstIndex()),
      leadingDimension_(rhs.leadingDimension())
{

}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorageView<T, Order, I, A>::~BandStorageView()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ElementType &
BandStorageView<T, Order, I, A>::operator()(IndexType row, IndexType col) const
{
    ASSERT(row>=firstIndex_);
    ASSERT(row<firstIndex_+numRows_);
    ASSERT(col>=firstIndex_);
    ASSERT(col<firstIndex_+numCols_);

    ASSERT(max(firstIndex_,col-numSuperDiags_) <= row);
    ASSERT(row <= min(numRows_+firstIndex_-1,col+numSubDiags_));

    if (Order == ColMajor) {
        const IndexType i = numSuperDiags_+row-col;
        const IndexType j = col-firstIndex_;
        return data_[j*leadingDimension_+i];
    }
    const IndexType i = numSubDiags_+col-row;
    const IndexType j = row-firstIndex_;
    return data_[j*leadingDimension_+i];
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ElementType &
BandStorageView<T, Order, I, A>::operator()(IndexType row, IndexType col)
{
    ASSERT(row>=firstIndex_);
    ASSERT(row<firstIndex_+numRows_);
    ASSERT(col>=firstIndex_);
    ASSERT(col<firstIndex_+numCols_);

    ASSERT(max(firstIndex_,col-numSuperDiags_) <= row);
    ASSERT(row <= min(numRows_+firstIndex_-1,col+numSubDiags_));
    if (Order == ColMajor) {
        const IndexType i = numSuperDiags_+row-col;
        const IndexType j = col-firstIndex_;
        return data_[j*leadingDimension_+i];
    }
    const IndexType i = numSubDiags_+col-row;
    const IndexType j = row-firstIndex_;
    return data_[j*leadingDimension_+i];
}

//-- methods -------------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::firstRow() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::lastRow() const
{
    return firstIndex_+numRows_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::firstCol() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::lastCol() const
{
    return firstIndex_+numCols_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::firstIndex() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::lastIndex() const
{
    return firstIndex_+numCols_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::numRows() const
{
    return numRows_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::numCols() const
{
    return numCols_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::dim() const
{
    ASSERT(numCols_ == numRows_);
    return numCols_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::numSubDiags() const
{
    return numSubDiags_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::numSuperDiags() const
{
    return numSuperDiags_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::leadingDimension() const
{
    return leadingDimension_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension()-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::IndexType
BandStorageView<T, Order, I, A>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()-1
                             : 1;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ElementType *
BandStorageView<T, Order, I, A>::data() const
{
    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ElementType *
BandStorageView<T, Order, I, A>::data()
{
    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::Allocator &
BandStorageView<T, Order, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorageView<T, Order, I, A>::resize(IndexType DEBUG_VAR(numRows),
                                        IndexType DEBUG_VAR(numCols),
                                        IndexType DEBUG_VAR(numSubDiags),
                                        IndexType DEBUG_VAR(numSuperDiags),
                                        IndexType firstIndex,
                                        const ElementType &)
{
    ASSERT(numRows_==numRows);
    ASSERT(numCols_==numCols);
    ASSERT(numSubDiags_==numSubDiags);
    ASSERT(numSuperDiags_==numSuperDiags);

    changeIndexBase(firstIndex);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename FS>
bool
BandStorageView<T, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.numRows(), rhs.numCols(),
                  rhs.numSubDiags(), rhs.numSuperDiags(), rhs.firstCol(),
                  value);
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorageView<T, Order, I, A>::fill(const ElementType &value)
{
    using std::max;
    using std::min;

    for (IndexType row = firstIndex_; row <= firstIndex_+numRows_-1; ++row)
    {
        const IndexType mincol = max(firstIndex_,row-numSubDiags_);
        const IndexType maxcol = min(row+numSuperDiags_,numCols_+firstIndex_-1);
        for (IndexType col = mincol; col <= maxcol; ++col)
            operator()(row, col) = value;
    }

    return true;
}


template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorageView<T, Order, I, A>::fillRandom()
{
    using std::max;
    using std::min;

    for (IndexType row = firstIndex_; row <= firstIndex_+numRows_-1; ++row)
    {
        const IndexType mincol = max(firstIndex_,row-numSubDiags_);
        const IndexType maxcol = min(row+numSuperDiags_,numCols_+firstIndex_-1);
        for (IndexType col = mincol; col <= maxcol; ++col)
            operator()(row, col) = randomValue<T>();
    }

    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorageView<T, Order, I, A>::changeIndexBase(IndexType firstIndex)
{
    firstIndex_ = firstIndex;
}

// view of fullstorage scheme as an array
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::arrayView(IndexType firstViewIndex) const
{

    if (Order==RowMajor) {
        return ConstArrayView(leadingDimension_*numRows_,
                              data_,
                              IndexType(1),
                              firstViewIndex,
                              allocator());
    }

    return ConstArrayView(leadingDimension_*numRows_,
                          data_,
                          IndexType(1),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::arrayView(IndexType firstViewIndex)
{
    if (Order==RowMajor) {
        return ArrayView(leadingDimension_*numCols_,
                          data_,
                          IndexType(1),
                          firstViewIndex,
                          allocator());
    }

    return ArrayView(leadingDimension_*numRows_,
                      data_,
                      IndexType(1),
                      firstViewIndex,
                      allocator());
}

// View a diagonal
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::viewDiag(IndexType diag,
                                          IndexType firstViewIndex) const
{

    ASSERT( diag <= numSuperDiags_);
    ASSERT(-diag <= numSubDiags_);

    using std::min;

    const IndexType i = (diag < 0) ? -diag+firstViewIndex: firstViewIndex;
    const IndexType j = (diag > 0) ?  diag+firstViewIndex: firstViewIndex;
    const IndexType length = (diag<=0) ? min(numCols_, numRows_+diag)
                                       : min(numCols_-diag, numRows_);

    return ConstArrayView(length-firstViewIndex+firstIndex_,
                          &(operator()(i, j)),
                          leadingDimension_,
                          firstIndex_, allocator_);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::viewDiag(IndexType diag,
                                          IndexType firstViewIndex)
{
    ASSERT( diag <= numSuperDiags_);
    ASSERT(-diag <= numSubDiags_);

    using std::min;

    const IndexType i = (diag < 0) ? -diag+firstViewIndex: firstViewIndex;
    const IndexType j = (diag > 0) ?  diag+firstViewIndex: firstViewIndex;
    const IndexType length = (diag<=0) ? min(numCols_, numRows_+diag)
                                       : min(numCols_-diag, numRows_);

    return ArrayView(length-firstViewIndex+firstIndex_,
                     &(operator()(i, j)),
                     leadingDimension_,
                     firstIndex_, allocator_);
}


template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstView
BandStorageView<T, Order, I, A>::viewDiags(IndexType fromDiag,
                                           IndexType toDiag) const
{
    ASSERT(fromDiag<=toDiag);
    IndexType numRows = numRows_;
    IndexType numCols = numCols_;

    if (fromDiag>0) {
        numCols = numCols_ - fromDiag;
        if (numRows_<numCols_) {
            if (numCols_-numRows_ < fromDiag) {
                numRows = numCols_ - fromDiag;
            }
        } else {
            numRows = numCols_ - fromDiag;
        }
    }
    if (toDiag<0) {
        numRows = numRows_ + toDiag;
        if (numCols_<numRows_) {
            if (numRows_-numCols_ < -toDiag) {
              numCols = numRows_ + toDiag;
            }
        } else {
            numCols = numRows_ + toDiag;
        }
    }

    const IndexType i = firstIndex_ - ((toDiag<0) ? toDiag : 0);
    const IndexType j = firstIndex_ + ((fromDiag>0) ? fromDiag : 0);

    if (Order == RowMajor ) {
        if (toDiag < 0) {
            return ConstView(numRows, numCols, -fromDiag+toDiag, 0,
                             &(operator()(i,j)) + fromDiag-toDiag,
                            leadingDimension_,
                             firstIndex_, allocator_);
        }
        if (fromDiag > 0) {
            return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                             &(operator()(i,j)),
                             leadingDimension_,
                             firstIndex_, allocator_);
        }
        return ConstView(numRows, numCols, -fromDiag, toDiag,
                         data_ + numSubDiags_+fromDiag,
                         leadingDimension_,
                         firstIndex_, allocator_);
    }

    if (toDiag < 0) {
        return ConstView(numRows, numCols, -fromDiag+toDiag, 0,
                         &(operator()(i,j)),
                         leadingDimension_,
                         firstIndex_, allocator_);
    }
    if (fromDiag > 0) {
        return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                         &(operator()(i,j)) + fromDiag-toDiag,
                         leadingDimension_,
                         firstIndex_, allocator_);
    }
    return ConstView(numRows_, numCols_, -fromDiag, toDiag,
                     data_ + (numSuperDiags_-toDiag),
                     leadingDimension_,
                     firstIndex_, allocator_);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::View
BandStorageView<T, Order, I, A>::viewDiags(IndexType fromDiag, IndexType toDiag)
{
    ASSERT(fromDiag<=toDiag);
    IndexType numRows = numRows_;
    IndexType numCols = numCols_;

    if (fromDiag>0) {
        numCols = numCols_ - fromDiag;
        if (numRows_<numCols_) {
            if (numCols_-numRows_ < fromDiag) {
                numRows = numCols_ - fromDiag;
            }
        } else {
            numRows = numCols_ - fromDiag;
        }
    }
    if (toDiag<0) {
        numRows = numRows_ + toDiag;
        if (numCols_<numRows_) {
            if (numRows_-numCols_ < -toDiag) {
              numCols = numRows_ + toDiag;
            }
        } else {
            numCols = numRows_ + toDiag;
        }
    }

    const IndexType i = firstIndex_ - ((toDiag<0) ? toDiag : 0);
    const IndexType j = firstIndex_ + ((fromDiag>0) ? fromDiag : 0);

    if (Order == RowMajor ) {
        if (toDiag < 0) {
            return View(numRows, numCols, -fromDiag+toDiag, 0,
                        &(operator()(i,j)) + fromDiag-toDiag ,
                        leadingDimension_,
                        firstIndex_, allocator_);
        }
        if (fromDiag > 0) {
            return View(numRows, numCols, 0, toDiag-fromDiag,
                        &(operator()(i,j)),
                        leadingDimension_,
                        firstIndex_, allocator_);
        }
        return View(numRows, numCols, -fromDiag, toDiag,
                    data_ + numSubDiags_+fromDiag,
                    leadingDimension_,
                    firstIndex_, allocator_);
    }

    if (toDiag < 0) {
        return View(numRows, numCols, -fromDiag+toDiag, 0,
                    &(operator()(i,j)),
                    leadingDimension_,
                    firstIndex_, allocator_);
    }
    if (fromDiag > 0) {
        return View(numRows, numCols, 0, toDiag-fromDiag,
                    &(operator()(i,j)) + fromDiag-toDiag,
                   leadingDimension_,
                    firstIndex_, allocator_);
    }
    return View(numRows_, numCols_, -fromDiag, toDiag,
                data_ + (numSuperDiags_-toDiag),
                leadingDimension_,
                firstIndex_, allocator_);
}


// view of single row
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::viewRow(IndexType row,
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

    IndexType length = 1
                     + min(row-1, numSubDiags())
                     + min(numRows()-row, numSuperDiags());

    return ConstArrayView(length,
                          &(operator()(row,
                                       firstCol()+max(0,row-1-numSubDiags()))),
                          strideCol(),
                          firstViewIndex,
                          allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::viewRow(IndexType row,
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

    IndexType length = 1
                     + min(row-1,numSubDiags())
                     + min(numRows()-row,numSuperDiags());

    return ArrayView(length,
                     &(operator()(row, firstCol()+max(0,row-1-numSubDiags()))),
                     strideCol(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::viewRow(IndexType row,
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
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::viewRow(IndexType row,
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

// view of single row
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::viewCol(IndexType col,
                                         IndexType firstViewIndex) const
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

    IndexType length = 1
                     + min(col-1, numSuperDiags())
                     + min(numCols()-col, numSubDiags());

    return ArrayView(length,
                     &(operator()(firstRow()+max(0,col-1-numSuperDiags()),
                                  col)),
                     strideRow(),
                     firstViewIndex,
                     allocator());

}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::viewCol(IndexType col,
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

    IndexType length = 1
                     + min(col-1, numSuperDiags())
                     + min(numCols()-col, numSubDiags());

    return ArrayView(length,
                     &(operator()(firstRow()+max(0,col-1-numSuperDiags()),
                                  col)),
                     strideRow(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorageView<T, Order, I, A>::ConstArrayView
BandStorageView<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
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
typename BandStorageView<T, Order, I, A>::ArrayView
BandStorageView<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
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
                     &(operator()(firstRow(), col)),
                     strideRow()*stride,
                     firstViewIndex,
                     allocator());
}



} // namespace flens

#endif // FLENS_STORAGE_BANDSTORAGE_BANDSTORAGEVIEW_TCC
