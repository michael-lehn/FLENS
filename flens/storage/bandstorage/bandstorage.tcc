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

namespace flens {

//= Constructors

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage()
    :  _data(0),
       _numRows(0), _numCols(0),
       _numSubDiags(0), _numSuperDiags(0),
       _firstIndex(I::defaultIndexBase)
{
}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage(IndexType numRows, IndexType numCols,
                                         IndexType numSubDiags, IndexType numSuperDiags,
                                         IndexType firstIndex,
                                         const ElementType &value,
                                         const Allocator &allocator)
    : _data(0), _allocator(allocator),
      _numRows(numRows), _numCols(numCols),
      _numSubDiags(numSubDiags), _numSuperDiags(numSuperDiags),
      _firstIndex(firstIndex)
{

    ASSERT(_numRows>=0);
    ASSERT(_numCols>=0);
    ASSERT(_numSubDiags>=0);
    ASSERT(_numSuperDiags>=0);

    _allocate(value);

}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::BandStorage(const BandStorage &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _numSubDiags(rhs.numSubDiags()), _numSuperDiags(rhs.numSuperDiags()),
      _firstIndex(rhs.firstIndex())
{
    _allocate(ElementType());
    const IndexType leadingDimension = _numSubDiags+_numSuperDiags+1;
    if (Order==ColMajor) {
        cxxblas::copy(leadingDimension*_numCols, rhs.data(), 1, _data, 1);
    }
    else {
        cxxblas::copy(leadingDimension*_numRows, rhs.data(), 1, _data, 1);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
template <typename RHS>
BandStorage<T, Order, I, A>::BandStorage(const RHS &rhs)
    : _data(0), _allocator(rhs.allocator()),
      _numRows(rhs.numRows()), _numCols(rhs.numCols()),
      _numSubDiags(rhs.numSubDiags()), _numSuperDiags(rhs.numSuperDiags()),
      _firstIndex(rhs.firstIndex())
{
    using std::max;
    using std::min;

    _allocate(ElementType());

    for (IndexType row = _firstIndex; row <= _firstIndex+_numRows-1; ++row)
    {
        const IndexType mincol = max(_firstIndex,row-_numSubDiags);
        const IndexType maxcol = min(row+_numSuperDiags,_numCols+_firstIndex-1);
        for (IndexType col = mincol; col <= maxcol; ++col)
            operator()(row, col) = rhs.operator()(row,col);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
BandStorage<T, Order, I, A>::~BandStorage()
{
    _release();
}

//-- operators -----------------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ElementType &
BandStorage<T, Order, I, A>::operator()(IndexType row, IndexType col) const
{

    ASSERT(row>=_firstIndex);
    ASSERT(row<_firstIndex+_numRows);
    ASSERT(col>=_firstIndex);
    ASSERT(col<_firstIndex+_numCols);

    ASSERT(max(_firstIndex,col-_numSuperDiags) <= row);
    ASSERT(row <= min(_numRows+_firstIndex-1,col+_numSubDiags));

    const IndexType leadingDimension = _numSubDiags+_numSuperDiags+1;
    if (Order == ColMajor) {
        const IndexType i = _numSuperDiags+row-col;
        const IndexType j = col-_firstIndex;

        return _data[j*leadingDimension+i];
    }

    const IndexType i = _numSubDiags+col-row;
    const IndexType j = row-_firstIndex;
    return _data[j*leadingDimension+i];
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ElementType &
BandStorage<T, Order, I, A>::operator()(IndexType row, IndexType col)
{
    using std::max;
    using std::min;
    ASSERT(row>=_firstIndex);
    ASSERT(row<_firstIndex+_numRows);
    ASSERT(col>=_firstIndex);
    ASSERT(col<_firstIndex+_numCols);

    ASSERT(max(_firstIndex,col-_numSuperDiags) <= row);
    ASSERT(row <= min(_numRows+_firstIndex-1,col+_numSubDiags));

    const IndexType leadingDimension = _numSubDiags+_numSuperDiags+1;
    if (Order == ColMajor) {
        const IndexType i = _numSuperDiags+row-col;
        const IndexType j = col-_firstIndex;
        return _data[j*leadingDimension+i];
    }

    const IndexType i = _numSubDiags+col-row;
    const IndexType j = row-_firstIndex;
    return _data[j*leadingDimension+i];

}

//-- Methods -------------------------------------------------------------------
template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstRow() const
{
    return _firstIndex;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastRow() const
{
    return _firstIndex+_numRows-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstCol() const
{
    return _firstIndex;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastCol() const
{
    return _firstIndex+_numCols-1;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::firstIndex() const
{
    return _firstIndex;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::lastIndex() const
{
    return _firstIndex+_numCols-1;
}


template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numRows() const
{
    return _numRows;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numCols() const
{
    return _numCols;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::dim() const
{
    ASSERT(_numCols == _numRows);
    return _numCols;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numSubDiags() const
{
    return _numSubDiags;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::numSuperDiags() const
{
    return _numSuperDiags;
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::IndexType
BandStorage<T, Order, I, A>::leadingDimension() const
{
    return std::max(_numSubDiags+_numSuperDiags+1, IndexType(1));
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

    return _data;
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

    return _data;
}



template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::Allocator &
BandStorage<T, Order, I, A>::allocator() const
{
    return _allocator;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorage<T, Order, I, A>::resize(IndexType numRows, IndexType numCols,
                                    IndexType numSubDiags, IndexType numSuperDiags,
                                    IndexType firstIndex,
                                    const ElementType &value)
{
    if ((_numSubDiags!=numSubDiags) ||(_numSuperDiags!=numSuperDiags)
      || (_numRows!=numRows) || (_numCols!=numCols)) {
        _release();
        _numSubDiags = numSubDiags,
        _numSuperDiags = numSuperDiags,
        _numRows = numRows;
        _numCols = numCols;
        _firstIndex = firstIndex;
        _allocate(value);
        return true;
    }
    _setIndexBase(firstIndex);
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
    const IndexType m = _numSubDiags+_numSuperDiags+1;
    if (Order==ColMajor) {
        std::fill_n(_data, m*_numCols, value);
    }
    else {
        std::fill_n(_data, m*_numRows, value);
    }

    return true;
}

template <typename T, StorageOrder Order, typename I, typename A>
bool
BandStorage<T, Order, I, A>::fillRandom()
{
    const IndexType m = _numSubDiags+_numSuperDiags+1;
    if (Order==ColMajor) {
        for (IndexType i=0; i<m*_numCols;++i) {
            _data[i] = randomValue<T>();
        }
    }
    else {
        for (IndexType i=0; i<m*_numRows;++i) {
            _data[i] = randomValue<T>();
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
        return ConstArrayView((_numSubDiags+_numSuperDiags+1)*_numRows,
                              _data,
                              IndexType(1),
                              firstViewIndex,
                              _allocator);
    }

    return ConstArrayView((_numSubDiags+_numSuperDiags+1)*_numRows,
                          _data,
                          IndexType(1),
                          firstViewIndex,
                          _allocator);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::arrayView(IndexType firstViewIndex)
{
    if (Order==RowMajor) {
        return ArrayView((_numSubDiags+_numSuperDiags+1)*_numRows,
                              _data,
                              IndexType(1),
                              firstViewIndex,
                              allocator());
    }

    return ArrayView((_numSubDiags+_numSuperDiags+1)*_numCols,
                          _data,
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

    ASSERT( diag <= _numSuperDiags);
    ASSERT(-diag <= _numSubDiags);

    using std::min;

    const IndexType i = (diag < 0) ? -diag+firstViewIndex: firstViewIndex;
    const IndexType j = (diag > 0) ?  diag+firstViewIndex: firstViewIndex;
    const IndexType length = (diag<=0) ? min(_numCols, _numRows+diag)
                                       : min(_numCols-diag, _numRows);

    return ConstArrayView(length-firstViewIndex+_firstIndex,
                          &(this->operator()(i, j)),
                          _numSubDiags+_numSuperDiags+1,
                          _firstIndex, _allocator);

}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::ArrayView
BandStorage<T, Order, I, A>::viewDiag(IndexType diag,
                                      IndexType firstViewIndex)
{
    ASSERT( diag <= _numSuperDiags);
    ASSERT(-diag <= _numSubDiags);

    using std::min;

    const IndexType i = (diag < 0) ? -diag+firstViewIndex: firstViewIndex;
    const IndexType j = (diag > 0) ?  diag+firstViewIndex: firstViewIndex;
    const IndexType length = (diag<=0) ? min(_numCols, _numRows+diag)
                                       : min(_numCols-diag, _numRows);

    return ArrayView(length-firstViewIndex+_firstIndex,
                     &(this->operator()(i, j)),
                     _numSubDiags+_numSuperDiags+1,
                     _firstIndex, _allocator);


}

// View some diagonals
template <typename T, StorageOrder Order, typename I, typename A>
const typename BandStorage<T, Order, I, A>::ConstView
BandStorage<T, Order, I, A>::viewDiags(IndexType fromDiag, IndexType toDiag) const
{

    ASSERT(fromDiag<=toDiag);
    IndexType numRows = _numRows;
    IndexType numCols = _numCols;

    if (fromDiag>0) {
        numCols = _numCols - fromDiag;
        if (_numRows<_numCols) {
            if (_numCols-_numRows < fromDiag) {
                numRows = _numCols - fromDiag;
            }
        } else {
            numRows = _numCols - fromDiag;
        }
    }
    if (toDiag<0) {
        numRows = _numRows + toDiag;
        if (_numCols<_numRows) {
            if (_numRows-_numCols < -toDiag) {
              numCols = _numRows + toDiag;
            }
        } else {
            numCols = _numRows + toDiag;
        }
    }

    const IndexType i = _firstIndex - ((toDiag<0) ? toDiag : 0);
    const IndexType j = _firstIndex + ((fromDiag>0) ? fromDiag : 0);

    if (Order == RowMajor ) {
        if (toDiag < 0) {
            return ConstView(numRows, numCols, -fromDiag+toDiag, 0,
                             &(this->operator()(i,j)) + fromDiag-toDiag ,
                             _numSubDiags+_numSuperDiags+1,
                             _firstIndex, _allocator);
        }
        if (fromDiag > 0) {
            return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                             &(this->operator()(i,j)),
                             _numSubDiags+_numSuperDiags+1,
                             _firstIndex, _allocator);
        }
        return ConstView(numRows, numCols, -fromDiag, toDiag,
                         _data + _numSubDiags+fromDiag,
                         _numSubDiags+_numSuperDiags+1,
                         _firstIndex, _allocator);
    }

    if (toDiag < 0) {
        return ConstView(numRows, numCols, -fromDiag+toDiag, 0,
                         &(this->operator()(i,j)),
                         _numSubDiags+_numSuperDiags+1,
                         _firstIndex, _allocator);
    }
    if (fromDiag > 0) {
        return ConstView(numRows, numCols, 0, toDiag-fromDiag,
                         &(this->operator()(i,j)) + fromDiag-toDiag,
                         _numSubDiags+_numSuperDiags+1,
                         _firstIndex, _allocator);
    }
    return ConstView(_numRows, _numCols, -fromDiag, toDiag,
                     _data + (_numSuperDiags-toDiag),
                     _numSubDiags+_numSuperDiags+1,
                     _firstIndex, _allocator);
}

template <typename T, StorageOrder Order, typename I, typename A>
typename BandStorage<T, Order, I, A>::View
BandStorage<T, Order, I, A>::viewDiags(IndexType fromDiag, IndexType toDiag)
{
    // WORKING !!!!
    ASSERT(fromDiag<=toDiag);
    IndexType numRows = _numRows;
    IndexType numCols = _numCols;

    if (fromDiag>0) {
        numCols = _numCols - fromDiag;
        if (_numRows<_numCols) {
            if (_numCols-_numRows < fromDiag) {
                numRows = _numCols - fromDiag;
            }
        } else {
            numRows = _numCols - fromDiag;
        }
    }
    if (toDiag<0) {
        numRows = _numRows + toDiag;
        if (_numCols<_numRows) {
            if (_numRows-_numCols < -toDiag) {
              numCols = _numRows + toDiag;
            }
        } else {
            numCols = _numRows + toDiag;
        }
    }

    const IndexType i = _firstIndex - ((toDiag<0) ? toDiag : 0);
    const IndexType j = _firstIndex + ((fromDiag>0) ? fromDiag : 0);

    if (Order == RowMajor ) {
        if (toDiag < 0) {
            return View(numRows, numCols, -fromDiag+toDiag, 0,
                        &(this->operator()(i,j)) + fromDiag-toDiag ,
                        _numSubDiags+_numSuperDiags+1,
                        _firstIndex, _allocator);
        }
        if (fromDiag > 0) {
            return View(numRows, numCols, 0, toDiag-fromDiag,
                        &(this->operator()(i,j)),
                        _numSubDiags+_numSuperDiags+1,
                        _firstIndex, _allocator);
        }
        return View(numRows, numCols, -fromDiag, toDiag,
                    _data + _numSubDiags+fromDiag,
                    _numSubDiags+_numSuperDiags+1,
                    _firstIndex, _allocator);
    }

    if (toDiag < 0) {
        return View(numRows, numCols, -fromDiag+toDiag, 0,
                    &(this->operator()(i,j)),
                    _numSubDiags+_numSuperDiags+1,
                    _firstIndex, _allocator);
    }
    if (fromDiag > 0) {
        return View(numRows, numCols, 0, toDiag-fromDiag,
                    &(this->operator()(i,j)) + fromDiag-toDiag,
                    _numSubDiags+_numSuperDiags+1,
                    _firstIndex, _allocator);
    }
    return View(_numRows, _numCols, -fromDiag, toDiag,
                _data + (_numSuperDiags-toDiag),
                _numSubDiags+_numSuperDiags+1,
                _firstIndex, _allocator);

}



//-- Private Methods -----------------------------------------------------------

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::_setIndexBase(IndexType firstIndex)
{
    _firstIndex = firstIndex;
}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::_raw_allocate()
{
    ASSERT(!_data);
    ASSERT(_numRows>0);
    ASSERT(_numCols>0);

    const IndexType m = _numSubDiags+_numSuperDiags+1;
    if (Order==ColMajor) {
        _data = _allocator.allocate(m*_numCols);
    }
    else {
        _data = _allocator.allocate(m*_numRows);
    }

    _setIndexBase(_firstIndex);

}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::_allocate(const ElementType &value)
{

    if (numRows()*numCols()==0) {
        return;
    }

    _raw_allocate();
    T *p = _data;

    IndexType numArrayElements = (_numSubDiags+_numSuperDiags+1)*_numCols;

     if (Order == RowMajor) {
        numArrayElements = (_numSubDiags+_numSuperDiags+1)*_numRows;
    }

    for (IndexType i=0; i<numArrayElements; ++i) {
        _allocator.construct(p++, value);
    }
}

template <typename T, StorageOrder Order, typename I, typename A>
void
BandStorage<T, Order, I, A>::_release()
{
    if (_data) {
        T *p = _data;
        IndexType numElements = (_numSubDiags+_numSuperDiags+1)*_numCols;
        if (Order == RowMajor)
            numElements = (_numSubDiags+_numSuperDiags+1)*_numRows;

        for (IndexType i=0; i<numElements; ++i) {
            _allocator.destroy(p++);
        }
         _allocator.deallocate(data(), numElements);
        _data = 0;
    }
    ASSERT(_data==0);
}

} // namespace flens

#endif // FLENS_STORAGE_BANDSTORAGE_BANDSTORAGE_TCC
