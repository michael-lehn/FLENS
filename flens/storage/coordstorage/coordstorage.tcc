/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_TCC
#define FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_TCC 1

#include <iostream>
#include <flens/storage/coordstorage/coordstorage.h>

namespace flens {

template <typename T, typename I, typename Cmp>
CoordStorage<T,I,Cmp>::CoordStorage(IndexType numRows, IndexType numCols,
                                    IndexType densityEstimate,
                                    IndexType firstRow,
                                    IndexType firstCol)
    : _numRows(numRows), _numCols(numCols),
      _firstRow(firstRow), _firstCol(firstCol),
      _lastSortedCoord(0), _isSorted(true), _isAccumulated(true)
{
    _coord.reserve(densityEstimate*_numRows);
}

template <typename T, typename I, typename Cmp>
CoordStorage<T,I,Cmp>::~CoordStorage()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I, typename Cmp>
typename CoordStorage<T,I,Cmp>::ElementProxy
CoordStorage<T,I,Cmp>::operator()(IndexType row, IndexType col)
{
    ASSERT(row>=_firstRow);
    ASSERT(row<_firstRow+_numRows);
    ASSERT(col>=_firstCol);
    ASSERT(col<_firstCol+_numCols);

    if (_coord.size()>=_coord.capacity()) {
        accumulate();
        _coord.reserve(_coord.capacity() + _numRows);
    }

    _coord.push_back(CoordType(row, col, ElementType(0)));
    _isAccumulated = false;

    size_t lastIndex = _coord.size()-1;
    if ((lastIndex>0) && _isSorted) {
        if (_less(_coord[lastIndex-1], _coord[lastIndex])) {
            _lastSortedCoord = lastIndex;
        } else {
            _isSorted = false;
        }
    }

    return &_coord[lastIndex].value;
}

//-- methods -------------------------------------------------------------------

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::firstRow() const
{
    return _firstRow;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::lastRow() const
{
    return _firstRow+_numRows-1;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::firstCol() const
{
    return _firstCol;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::lastCol() const
{
    return _firstCol+_numCols-1;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::numRows() const
{
    return _numRows;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::numCols() const
{
    return _numCols;
}

template <typename T, typename I, typename Cmp>
void
CoordStorage<T,I,Cmp>::accumulate() const
{
//
//  Quick return if possible
//
    if (_isAccumulated) {
        ASSERT(_isSorted);
        return;
    }

//
//  sort
//
    if (!_isSorted) {
        std::sort(_coord.begin()+_lastSortedCoord+1, _coord.end(), _less);
    }
    if ((_lastSortedCoord>0) && (_lastSortedCoord<_coord.size())) {
        std::inplace_merge(_coord.begin(),
                           _coord.begin() + _lastSortedCoord+1,
                           _coord.end(),
                           _less);
    }
    _isSorted = true;

//
//  accumulate values
//
    size_t k, K;
    for (k=0, K=1; K<_coord.size(); ++k, ++K) {
        while ((K<_coord.size()) && (!_less(_coord[k], _coord[K]))) {
            _coord[k].value += _coord[K].value;
            _coord[K].value = ElementType(0);
            ++K;
        }
        if (K<_coord.size()) {
            _coord[k+1] = _coord[K];
        }
    }
    if ((k<_coord.size()) && (K-k-1>0)) {
#       ifndef NDEBUG
        size_t oldCapacity = _coord.capacity();
#       endif

        _coord.erase(_coord.end()-(K-k-1), _coord.end());

#       ifndef NDEBUG
        if (oldCapacity!=_coord.capacity()) {
            std::cerr << "[WARNING] Possible performance bottleneck in "
                      << "CoordStorage<T,I,Cmp>::accumulate()"
                      << std::endl;
        }
#       endif
    }
    _lastSortedCoord = _coord.size()-1;
    _isAccumulated = true;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::CoordVector &
CoordStorage<T,I,Cmp>::coordVector() const
{
    return _coord;
}

template <typename T, typename I, typename Cmp>
const typename CoordStorage<T,I,Cmp>::IndexType
CoordStorage<T,I,Cmp>::numNonZeros() const
{
    return _coord.size();
}


//-- Coord ---------------------------------------------------------------------

template <typename T, typename IndexType>
Coord<T,IndexType>::Coord(IndexType _row, IndexType _col, const T &_value)
    : row(_row), col(_col), value(_value)
{
}

//-- CoordProxy ----------------------------------------------------------------

template <typename T>
CoordProxy<T>::CoordProxy(T *_value)
    : value(_value)
{
}

template <typename T>
void
CoordProxy<T>::operator+=(const T &x)
{
    *value = x;
}

template <typename T>
void
CoordProxy<T>::operator-=(const T &x)
{
    *value = -x;
}


//-- CoordRowColCmp ------------------------------------------------------------

template <typename T, typename IndexType>
bool
CoordRowColCmp::operator()(const Coord<T, IndexType> &a,
                           const Coord<T, IndexType> &b) const
{
    if (a.row<b.row) {
        return true;
    }
    if ((a.row==b.row) && (a.col<b.col)) {
        return true;
    }
    return false;
}

//-- CoordColRowCmp ------------------------------------------------------------

template <typename T, typename IndexType>
bool
CoordColRowCmp::operator()(const Coord<T, IndexType> &a,
                           const Coord<T, IndexType> &b) const
{
    if (a.col<b.col) {
        return true;
    }
    if ((a.col==b.col) && (a.row<b.row)) {
        return true;
    }
    return false;
}


} // namespace flens

#endif // FLENS_STORAGE_COORDSTORAGE_COORDSTORAGE_TCC
