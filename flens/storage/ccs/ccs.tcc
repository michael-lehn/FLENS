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

#include <vector>
#include <flens/storage/ccs/ccs.h>

namespace flens {

template <typename T, typename I>
CCS<T,I>::CCS()
    : _numRows(0), _numCols(0),
      _indexBase(I::defaultIndexBase)
{
}

template <typename T, typename I>
CCS<T,I>::~CCS()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I>
void
CCS<T,I>::operator=(const CoordinateStorage &coordinateStorage)
{
    _compress(coordinateStorage);
}

//-- methods -------------------------------------------------------------------

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::indexBase() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::firstRow() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::lastRow() const
{
    return _indexBase+_numRows-1;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::firstCol() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::lastCol() const
{
    return _indexBase+_numCols-1;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numRows() const
{
    return _numRows;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numCols() const
{
    return _numCols;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexType
CCS<T,I>::numNonZeros() const
{
    return _values.length();
}

template <typename T, typename I>
const typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::rows() const
{
    return _rows;
}

template <typename T, typename I>
typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::rows()
{
    return _rows;
}

template <typename T, typename I>
const typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::cols() const
{
    return _cols;
}

template <typename T, typename I>
typename CCS<T,I>::IndexTypeVector &
CCS<T,I>::cols()
{
    return _cols;
}

template <typename T, typename I>
const typename CCS<T,I>::ElementTypeVector &
CCS<T,I>::values() const
{
    return _values;
}

template <typename T, typename I>
typename CCS<T,I>::ElementTypeVector &
CCS<T,I>::values()
{
    return _values;
}

template <typename T, typename I>
void
CCS<T,I>::_compress(const CoordinateStorage &coordinateStorage)
{
    _numRows  = coordinateStorage.numRows();
    _numCols  = coordinateStorage.numCols();
    _indexBase = coordinateStorage.indexBase();

//
//  Accumulate coords and get number of non zeros
//
    coordinateStorage.accumulate();
    IndexType nnz = coordinateStorage.numNonZeros();

//
//  Allocate memory for the CCS format
//
    _cols.resize(_numCols+1, _indexBase);
    _rows.resize(nnz, _indexBase);
    _values.resize(nnz, _indexBase);

    const auto &coord = coordinateStorage.coordVector();

    IndexType c = _indexBase;
    _cols(c) = _indexBase;

    for (size_t k=0; k<coord.size(); ++k) {
        while (coord[k].col>c) {
            _cols(c+1) = _indexBase + k;
            ++c;
        }
        _rows(_indexBase+k)   = coord[k].row;
        _values(_indexBase+k) = coord[k].value;
    }
    while (c<=lastCol()) {
        _cols(c+1) = _indexBase + coord.size();
        ++c;
    }
}

} // namespace flens

#endif // FLENS_STORAGE_CCS_CCS_TCC
