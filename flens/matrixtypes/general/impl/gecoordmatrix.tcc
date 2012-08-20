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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/general/impl/gecoordmatrix.h>

namespace flens {

template <typename CS>
GeCoordMatrix<CS>::GeCoordMatrix(IndexType numRows, IndexType numCols)
    : _engine(numRows, numCols)
{
}

template <typename CS>
GeCoordMatrix<CS>::GeCoordMatrix(IndexType numRows, IndexType numCols,
                                 IndexType firstRow, IndexType firstCol)
    : _engine(numRows, numCols, firstRow, firstCol)
{
}

template <typename CS>
GeCoordMatrix<CS>::GeCoordMatrix(const Range<IndexType> &rowRange,
                                 const Range<IndexType> &colRange)
    : _engine(rowRange.numTicks(), colRange.numTicks(),
              rowRange.firstIndex(), colRange.firstIndex())
{
    ASSERT(rowRange.stride()==1);
    ASSERT(colRange.stride()==1);
}

// -- operators ----------------------------------------------------------------
template <typename CS>
typename GeCoordMatrix<CS>::ElementProxy
GeCoordMatrix<CS>::operator()(IndexType row, IndexType col)
{
    return _engine(row, col);
}

// -- methods ------------------------------------------------------------------
template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::numRows() const
{
    return _engine.numRows();
}

template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::numCols() const
{
    return _engine.numCols();
}

template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename CS>
const typename GeCoordMatrix<CS>::IndexType
GeCoordMatrix<CS>::lastCol() const
{
    return _engine.lastCol();
}

// -- implementation -----------------------------------------------------------
template <typename CS>
const typename GeCoordMatrix<CS>::Engine &
GeCoordMatrix<CS>::engine() const
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_TCC
