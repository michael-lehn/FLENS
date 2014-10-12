/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_STORAGE_CCS_CCS_TCC
#define FLENS_STORAGE_CCS_CCS_TCC 1

#include <flens/storage/ccs/ccs.h>

namespace flens {

template <typename T, typename I>
CCS<T,I>::CCS()
    : numRows_(0), numCols_(0),
      indexBase_(I::defaultIndexBase)
{
}

template <typename T, typename I>
CCS<T,I>::~CCS()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I>
template <typename T2, typename I2>
void
CCS<T,I>::operator=(const CoordStorage<T2, CoordColRowCmp, I2> &coordStorage)
{
    compress_(coordStorage);
}

//-- methods -------------------------------------------------------------------

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::indexBase() const
{
    return indexBase_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::firstRow() const
{
    return indexBase_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::lastRow() const
{
    return indexBase_+numRows_-1;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::firstCol() const
{
    return indexBase_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::lastCol() const
{
    return indexBase_+numCols_-1;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numRows() const
{
    return numRows_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numCols() const
{
    return numCols_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numNonZeros() const
{
    return values_.length();
}

template <typename T, typename I>
const typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::rows() const
{
    return rows_;
}

template <typename T, typename I>
typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::rows()
{
    return rows_;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::cols() const
{
    return cols_;
}

template <typename T, typename I>
typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::cols()
{
    return cols_;
}

template <typename T, typename I>
const typename CCS<T,I>::ElementTypeVector &
CCS<T,I>::values() const
{
    return values_;
}

template <typename T, typename I>
typename CCS<T,I>::ElementTypeVector &
CCS<T,I>::values()
{
    return values_;
}

template <typename T, typename I>
template <typename T2, typename I2>
void
CCS<T,I>::compress_(const CoordStorage<T2, CoordColRowCmp, I2> &coordStorage)
{
    numRows_  = coordStorage.numRows();
    numCols_  = coordStorage.numCols();
    indexBase_ = coordStorage.indexBase();

//
//  Accumulate coords and get number of non zeros
//
    coordStorage.accumulate();
    IndexType nnz = coordStorage.numNonZeros();

//
//  Allocate memory for the CCS format
//
    cols_.resize(numCols_+1, indexBase_);
    rows_.resize(nnz, indexBase_);
    values_.resize(nnz, indexBase_);

    const auto &coord = coordStorage.coordVector();

    IndexType c = indexBase_;
    cols_(c) = indexBase_;

    for (size_t k=0; k<coord.size(); ++k) {
        while (coord[k].col>c) {
            cols_(c+1) = indexBase_ + k;
            ++c;
        }
        rows_(indexBase_+k)   = coord[k].row;
        values_(indexBase_+k) = coord[k].value;
    }
    while (c<=lastCol()) {
        cols_(c+1) = indexBase_ + coord.size();
        ++c;
    }
}

} // namespace flens

#endif // FLENS_STORAGE_CCS_CCS_TCC
