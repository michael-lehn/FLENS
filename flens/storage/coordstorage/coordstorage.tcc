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

#include <cxxstd/iostream.h>
#include <flens/storage/coordstorage/coordstorage.h>

namespace flens {

template <typename T, typename Cmp, typename I>
CoordStorage<T,Cmp,I>::CoordStorage(IndexType numRows, IndexType numCols,
                                    IndexType densityEstimate,
                                    IndexType indexBase)
    : numRows_(numRows), numCols_(numCols), indexBase_(indexBase),
      lastSortedCoord_(0), isSorted_(true), isAccumulated_(true)
{
    coord_.reserve(densityEstimate*numRows_);
}

template <typename T, typename Cmp, typename I>
CoordStorage<T,Cmp,I>::~CoordStorage()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename Cmp, typename I>
typename CoordStorage<T,Cmp,I>::ElementProxy
CoordStorage<T,Cmp,I>::operator()(IndexType row, IndexType col)
{
    ASSERT(row>=indexBase_);
    ASSERT(row<indexBase_+numRows_);
    ASSERT(col>=indexBase_);
    ASSERT(col<indexBase_+numCols_);

    if (coord_.size()>=coord_.capacity()) {
        accumulate();
        coord_.reserve(coord_.capacity() + numRows_);
    }

    coord_.push_back(CoordType(row, col, ElementType(0)));
    isAccumulated_ = false;

    size_t lastIndex = coord_.size()-1;
    if ((lastIndex>0) && isSorted_) {
        if (less_(coord_[lastIndex-1], coord_[lastIndex])) {
            lastSortedCoord_ = lastIndex;
        } else {
            isSorted_ = false;
        }
    }

    return &coord_[lastIndex].value;
}

//-- methods -------------------------------------------------------------------

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::indexBase() const
{
    return indexBase_;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::firstRow() const
{
    return indexBase_;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::lastRow() const
{
    return indexBase_+numRows_-1;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::firstCol() const
{
    return indexBase_;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::lastCol() const
{
    return indexBase_+numCols_-1;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::numRows() const
{
    return numRows_;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::numCols() const
{
    return numCols_;
}

template <typename T, typename Cmp, typename I>
void
CoordStorage<T,Cmp,I>::accumulate() const
{
//
//  Quick return if possible
//
    if (isAccumulated_) {
        ASSERT(isSorted_);
        return;
    }

//
//  sort
//
    if (!isSorted_) {
        std::sort(coord_.begin()+lastSortedCoord_+1, coord_.end(), less_);
    }
    if (lastSortedCoord_<coord_.size()) {
        std::inplace_merge(coord_.begin(),
                           coord_.begin() + lastSortedCoord_+1,
                           coord_.end(),
                           less_);
    }
    isSorted_ = true;

//
//  accumulate values
//
    size_t k, K;
    for (k=0, K=1; K<coord_.size(); ++k, ++K) {
        while ((K<coord_.size()) && (!less_(coord_[k], coord_[K]))) {
            coord_[k].value += coord_[K].value;
            coord_[K].value = ElementType(0);
            ++K;
        }
        if (K<coord_.size()) {
            coord_[k+1] = coord_[K];
        }
    }
    if ((k<coord_.size()) && (K-k-1>0)) {
#       ifndef NDEBUG
        size_t oldCapacity = coord_.capacity();
#       endif

        coord_.erase(coord_.end()-(K-k-1), coord_.end());

#       ifndef NDEBUG
        if (oldCapacity!=coord_.capacity()) {
            std::cerr << "[WARNING] Possible performance bottleneck in "
                      << "CoordStorage<T,Cmp,I>::accumulate()"
                      << std::endl;
        }
#       endif
    }
    lastSortedCoord_ = coord_.size()-1;
    isAccumulated_ = true;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::CoordVector &
CoordStorage<T,Cmp,I>::coordVector() const
{
    return coord_;
}

template <typename T, typename Cmp, typename I>
const typename CoordStorage<T,Cmp,I>::IndexType
CoordStorage<T,Cmp,I>::numNonZeros() const
{
    return coord_.size();
}


//-- Coord ---------------------------------------------------------------------

template <typename T, typename IndexType>
Coord<T,IndexType>::Coord(IndexType row_, IndexType col_, const T &value_)
    : row(row_), col(col_), value(value_)
{
}

//-- CoordProxy ----------------------------------------------------------------

template <typename T>
CoordProxy<T>::CoordProxy(T *value_)
    : value(value_)
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
