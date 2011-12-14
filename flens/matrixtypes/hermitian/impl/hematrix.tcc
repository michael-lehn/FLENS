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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_TCC
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_TCC 1

#include <flens/blas/level1/copy.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
HeMatrix<FS>::HeMatrix()
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(IndexType dim)
    : _engine(dim, dim)
{
    ASSERT(dim>=0);
}

template <typename FS>
HeMatrix<FS>::HeMatrix(IndexType dim, IndexType firstRow, IndexType firstCol)
    : _engine(dim, dim, firstRow, firstCol)
{
    ASSERT(dim>=0);
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const Engine &engine, StorageUpLo upLo)
    : _engine(engine), _upLo(upLo)
{
    ASSERT(_engine.numRows()==_engine.numCols());
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const HeMatrix &rhs)
    : HermitianMatrix<HeMatrix<FS> >(),
      _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const HeMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(HeMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
const typename HeMatrix<FS>::ElementType &
HeMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (_upLo==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
#   endif
    return _engine(row, col);
}

template <typename FS>
typename HeMatrix<FS>::ElementType &
HeMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (_upLo==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
#   endif
    return _engine(row, col);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::general() const
{
    return _engine;
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::general()
{
    return _engine;
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::dim() const
{
    ASSERT(_engine.numRows()==_engine.numCols());

    return _engine.numRows();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::lastCol() const
{
    return _engine.lastCol();
}

template <typename FS>
const typename HeMatrix<FS>::ElementType *
HeMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename HeMatrix<FS>::ElementType *
HeMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename HeMatrix<FS>::Engine &
HeMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename HeMatrix<FS>::Engine &
HeMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
HeMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
StorageUpLo &
HeMatrix<FS>::upLo()
{
    return _upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_TCC
