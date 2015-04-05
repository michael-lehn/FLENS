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

#ifndef FLENS_STORAGE_BANDSTORAGE_BANDSTORAGE_TCC
#define FLENS_STORAGE_BANDSTORAGE_BANDSTORAGE_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/storage/bandstorage/bandstorage.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

namespace flens {

//= Constructors

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage()
    :  data_(0),
       numRows_(0), numCols_(0),
       numSubDiags_(0), numSuperDiags_(0),
       firstIndex_(I::defaultIndexBase)
{
}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage(IndexType numRows, IndexType numCols,
                                         IndexType numSubDiags,
                                         IndexType numSuperDiags,
                                         IndexType firstIndex,
                                         const ElementType &value,
                                         const Allocator &allocator)
    : data_(0), allocator_(allocator),
      numRows_(numRows), numCols_(numCols),
      numSubDiags_(numSubDiags), numSuperDiags_(numSuperDiags),
      firstIndex_(firstIndex)
{

    ASSERT(numRows_>=0);
    ASSERT(numCols_>=0);
    ASSERT(numSubDiags_>=0);
    ASSERT(numSuperDiags_>=0);

    allocate_(value);

}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage(const BandStorage &rhs)
    : data_(0), allocator_(rhs.allocator()),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      numSubDiags_(rhs.numSubDiags()), numSuperDiags_(rhs.numSuperDiags()),
      firstIndex_(rhs.firstIndex())
{
    allocate_(ElementType());
    const IndexType leadingDimension = numSubDiags_+numSuperDiags_+1;
    if (Order==ColMajor) {
        cxxblas::copy(leadingDimension*numCols_, rhs.data(), 1, data_, 1);
    }
    else {
        cxxblas::copy(leadingDimension*numRows_, rhs.data(), 1, data_, 1);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
BandStorage<T, Order, I, A>::BandStorage(const RHS &rhs)
    : data_(0), allocator_(rhs.allocator()),
      numRows_(rhs.numRows()), numCols_(rhs.numCols()),
      numSubDiags_(rhs.numSubDiags()), numSuperDiags_(rhs.numSuperDiags()),
      firstIndex_(rhs.firstIndex())
{
    using std::max;
    using std::min;

    allocate_(ElementType());

    for (IndexType row = firstIndex_; row <= firstIndex_+numRows_-1; ++row)
    {
        const IndexType mincol = max(firstIndex_,row-numSubDiags_);
        const IndexType maxcol = min(row+numSuperDiags_,numCols_+firstIndex_-1);
        for (IndexType col = mincol; col <= maxcol; ++col)
            operator()(row, col) = rhs.operator()(row,col);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::~BandStorage()
{
    release_();
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ElementType &
BandStorage<T, Order, I, A>::operator()(IndexType row, IndexType col) const
{

    ASSERT(row>=firstIndex_);
    ASSERT(row<firstIndex_+numRows_);
    ASSERT(col>=firstIndex_);
    ASSERT(col<firstIndex_+numCols_);

    ASSERT(max(firstIndex_,col-numSuperDiags_) <= row);
    ASSERT(row <= min(numRows_+firstIndex_-1,col+numSubDiags_));

    const IndexType leadingDimension = numSubDiags_+numSuperDiags_+1;
    if (Order == ColMajor) {
        const IndexType i = numSuperDiags_+row-col;
        const IndexType j = col-firstIndex_;

        return data_[j*leadingDimension+i];
    }

    const IndexType i = numSubDiags_+col-row;
    const IndexType j = row-firstIndex_;
    return data_[j*leadingDimension+i];
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ElementType &
BandStorage<T, Order, I, A>::operator()(IndexType row, IndexType col)
{
    using std::max;
    using std::min;
    ASSERT(row>=firstIndex_);
    ASSERT(row<firstIndex_+numRows_);
    ASSERT(col>=firstIndex_);
    ASSERT(col<firstIndex_+numCols_);

    ASSERT(max(firstIndex_,col-numSuperDiags_) <= row);
    ASSERT(row <= min(numRows_+firstIndex_-1,col+numSubDiags_));

    const IndexType leadingDimension = numSubDiags_+numSuperDiags_+1;
    if (Order == ColMajor) {
        const IndexType i = numSuperDiags_+row-col;
        const IndexType j = col-firstIndex_;
        return data_[j*leadingDimension+i];
    }

    const IndexType i = numSubDiags_+col-row;
    const IndexType j = row-firstIndex_;
    return data_[j*leadingDimension+i];

}

//-- Methods -------------------------------------------------------------------
template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstRow() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastRow() const
{
    return firstIndex_+numRows_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstCol() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastCol() const
{
    return firstIndex_+numCols_-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstIndex() const
{
    return firstIndex_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastIndex() const
{
    return firstIndex_+numCols_-1;
}


template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numRows() const
{
    return numRows_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numCols() const
{
    return numCols_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::dim() const
{
    ASSERT(numCols_ == numRows_);
    return numCols_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numSubDiags() const
{
    return numSubDiags_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numSuperDiags() const
{
    return numSuperDiags_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::leadingDimension() const
{
    return std::max(numSubDiags_+numSuperDiags_+1, IndexType(1));
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::strideRow() const
{
    return (Order==ColMajor) ? 1
                             : leadingDimension()-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::strideCol() const
{
    return (Order==ColMajor) ? leadingDimension()-1
                             : 1;
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ElementType *
BandStorage<T, Order, I, A>::data() const
{
#   ifndef NDEBUG
    if ((numRows()==0) || numCols()==0) {
        return 0;
    }
#   endif

    return data_;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ElementType *
BandStorage<T, Order, I, A>::data()
{
#   ifndef NDEBUG
    if ((numRows()==0) || numCols()==0) {
        return 0;
    }
#   endif

    return data_;
}



template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::Allocator &
BandStorage<T, Order, I, A>::allocator() const
{
    return allocator_;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorage<T, Order, I, A>::resize(IndexType numRows, IndexType numCols,
                                    IndexType numSubDiags,
                                    IndexType numSuperDiags,
                                    IndexType firstIndex,
                                    const ElementType &value)
{
    if ((numSubDiags_!=numSubDiags) ||(numSuperDiags_!=numSuperDiags)
      || (numRows_!=numRows) || (numCols_!=numCols)) {
        release_();
        numSubDiags_ = numSubDiags,
        numSuperDiags_ = numSuperDiags,
        numRows_ = numRows;
        numCols_ = numCols;
        firstIndex_ = firstIndex;
        allocate_(value);
        return true;
    }
    setIndexBase_(firstIndex);
    return false;
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename FS>
bool
BandStorage<T, Order, I, A>::resize(const FS &rhs, const ElementType &value)
{
    return resize(rhs.numRows(), rhs.numCols(),
                  rhs.numSubDiags(), rhs.numSuperDiags(),
                  rhs.firstIndex(),
                  value);
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorage<T, Order, I, A>::fill(const ElementType &value)
{
    const IndexType m = numSubDiags_+numSuperDiags_+1;
    if (Order==ColMajor) {
        std::fill_n(data_, m*numCols_, value);
    }
    else {
        std::fill_n(data_, m*numRows_, value);
    }

    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorage<T, Order, I, A>::fillRandom()
{
    const IndexType m = numSubDiags_+numSuperDiags_+1;
    if (Order==ColMajor) {
        for (IndexType i=0; i<m*numCols_;++i) {
            data_[i] = randomValue<T>();
        }
    }
    else {
        for (IndexType i=0; i<m*numRows_;++i) {
            data_[i] = randomValue<T>();
        }
    }

    return true;
}

// view of fullstorage scheme as an array
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::arrayView(IndexType firstViewIndex) const
{
    if (Order==RowMajor) {
        return ConstArrayView((numSubDiags_+numSuperDiags_+1)*numRows_,
                              data_,
                              IndexType(1),
                              firstViewIndex,
                              allocator_);
    }

    return ConstArrayView((numSubDiags_+numSuperDiags_+1)*numRows_,
                          data_,
                          IndexType(1),
                          firstViewIndex,
                          allocator_);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::arrayView(IndexType firstViewIndex)
{
    if (Order==RowMajor) {
        return ArrayView((numSubDiags_+numSuperDiags_+1)*numRows_,
                              data_,
                              IndexType(1),
                              firstViewIndex,
                              allocator());
    }

    return ArrayView((numSubDiags_+numSuperDiags_+1)*numCols_,
                          data_,
                          IndexType(1),
                          firstViewIndex,
                          allocator());
}


// view of a diagonal
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::viewDiag(IndexType diag,
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
                          numSubDiags_+numSuperDiags_+1,
                          firstIndex_, allocator_);

}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewDiag(IndexType diag,
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
                     numSubDiags_+numSuperDiags_+1,
                     firstIndex_, allocator_);


}

// View some diagonals
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstView
BandStorage<T, Order, I, A>::viewDiags(IndexType fromDiag,
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
                             numSubDiags_+numSuperDiags_+1,
                             firstIndex_, allocator_);
        }
        if (fromDiag > 0) {
            return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                             &(operator()(i,j)),
                             numSubDiags_+numSuperDiags_+1,
                             firstIndex_, allocator_);
        }
        return ConstView(numRows, numCols, -fromDiag, toDiag,
                         data_ + numSubDiags_+fromDiag,
                         numSubDiags_+numSuperDiags_+1,
                         firstIndex_, allocator_);
    }

    if (toDiag < 0) {
        return ConstView(numRows, numCols, -fromDiag+toDiag, 0,
                         &(operator()(i,j)),
                         numSubDiags_+numSuperDiags_+1,
                         firstIndex_, allocator_);
    }
    if (fromDiag > 0) {
        return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                         &(operator()(i,j)) + fromDiag-toDiag,
                         numSubDiags_+numSuperDiags_+1,
                         firstIndex_, allocator_);
    }
    return ConstView(numRows_, numCols_, -fromDiag, toDiag,
                     data_ + (numSuperDiags_-toDiag),
                     numSubDiags_+numSuperDiags_+1,
                     firstIndex_, allocator_);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::View
BandStorage<T, Order, I, A>::viewDiags(IndexType fromDiag, IndexType toDiag)
{
    // WORKING !!!!
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
                        numSubDiags_+numSuperDiags_+1,
                        firstIndex_, allocator_);
        }
        if (fromDiag > 0) {
            return View(numRows, numCols, 0, toDiag-fromDiag,
                        &(operator()(i,j)),
                        numSubDiags_+numSuperDiags_+1,
                        firstIndex_, allocator_);
        }
        return View(numRows, numCols, -fromDiag, toDiag,
                    data_ + numSubDiags_+fromDiag,
                    numSubDiags_+numSuperDiags_+1,
                    firstIndex_, allocator_);
    }

    if (toDiag < 0) {
        return View(numRows, numCols, -fromDiag+toDiag, 0,
                    &(operator()(i,j)),
                    numSubDiags_+numSuperDiags_+1,
                    firstIndex_, allocator_);
    }
    if (fromDiag > 0) {
        return View(numRows, numCols, 0, toDiag-fromDiag,
                    &(operator()(i,j)) + fromDiag-toDiag,
                    numSubDiags_+numSuperDiags_+1,
                    firstIndex_, allocator_);
    }
    return View(numRows_, numCols_, -fromDiag, toDiag,
                data_ + (numSuperDiags_-toDiag),
                numSubDiags_+numSuperDiags_+1,
                firstIndex_, allocator_);
}

// view of single row
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::viewRow(IndexType row,
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
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewRow(IndexType row,
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
                     + min(row-1, numSubDiags())
                     + min(numRows()-row, numSuperDiags());

    return ArrayView(length,
                     &(operator()(row, firstCol()+max(0,row-1-numSubDiags()))),
                     strideCol(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::viewRow(IndexType row,
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
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewRow(IndexType row,
                                     IndexType firstCol, IndexType lastCol,
                                     IndexType stride, IndexType firstViewIndex)
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
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::viewCol(IndexType col,
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
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewCol(IndexType col,
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
                     + min(col-1,numSuperDiags())
                     + min(numCols()-col,numSubDiags());

    return ArrayView(length,
                     &(operator()(firstRow()+max(0,col-1-numSuperDiags()),
                                  col)),
                     strideRow(),
                     firstViewIndex,
                     allocator());
}

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstArrayView
BandStorage<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
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
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewCol(IndexType firstRow, IndexType lastRow,
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

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstFullStorageView
BandStorage<T, Order, I, A>::viewFullStorage() const
{
    return ConstFullStorageView(numSubDiags()+numSuperDiags()+1,
                                max(numRows(),numCols()),
                                leadingDimension(),
                                data());
}


template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::FullStorageView
BandStorage<T, Order, I, A>::viewFullStorage()
{
    return FullStorageView(numSubDiags()+numSuperDiags()+1,
                           max(numRows(),numCols()),
                           leadingDimension(),
                           data());
}

//-- Private Methods -----------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::setIndexBase_(IndexType firstIndex)
{
    firstIndex_ = firstIndex;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::raw_allocate_()
{
    ASSERT(!data_);
    ASSERT(numRows_>0);
    ASSERT(numCols_>0);

    const IndexType m = numSubDiags_+numSuperDiags_+1;
    if (Order==ColMajor) {
        data_ = allocator_.allocate(m*numCols_);
    }
    else {
        data_ = allocator_.allocate(m*numRows_);
    }

    setIndexBase_(firstIndex_);

}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::allocate_(const ElementType &value)
{

    if (numRows()*numCols()==0) {
        return;
    }

    raw_allocate_();
    T *p = data_;

    IndexType numArrayElements = (numSubDiags_+numSuperDiags_+1)*numCols_;

     if (Order == RowMajor) {
        numArrayElements = (numSubDiags_+numSuperDiags_+1)*numRows_;
    }

    for (IndexType i=0; i<numArrayElements; ++i) {
        allocator_.construct(p++, value);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::release_()
{
    if (data_) {
        T *p = data_;
        IndexType numElements = (numSubDiags_+numSuperDiags_+1)*numCols_;
        if (Order == RowMajor)
            numElements = (numSubDiags_+numSuperDiags_+1)*numRows_;

        for (IndexType i=0; i<numElements; ++i) {
            allocator_.destroy(p++);
        }
         allocator_.deallocate(data(), numElements);
        data_ = 0;
    }
    ASSERT(data_==0);
}

} // namespace flens

#endif // FLENS_STORAGE_BANDSTORAGE_BANDSTORAGE_TCC
