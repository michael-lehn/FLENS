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
void
HeMatrix<FS>::operator=(const ElementType &value)
{
    engine().fill(_upLo, value);
}

template <typename FS>
HeMatrix<FS> &
HeMatrix<FS>::operator=(const HeMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
HeMatrix<FS> &
HeMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

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

// rectangular views
template <typename FS>
const typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols) const
{
    return general()(rows, cols);
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols)
{
    return general()(rows, cols);
}

// rectangular views (all rows selected)
template <typename FS>
const typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols) const
{
    return general()(_,cols);
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols)
{
    return general()(_,cols);
}

// rectangular views (all columns selected)
template <typename FS>
const typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_) const
{
    return general()(rows,_);
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_)
{
    return general()(rows,_);
}

// row view (vector view)
template <typename FS>
const typename HeMatrix<FS>::ConstVectorView
HeMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_) const
{
    return general()(row,_);
}

template <typename FS>
typename HeMatrix<FS>::VectorView
HeMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_)
{
    return general()(row,_);
}

template <typename FS>
const typename HeMatrix<FS>::ConstVectorView
HeMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols) const
{
    return general()(row,cols);
}

template <typename FS>
typename HeMatrix<FS>::VectorView
HeMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols)
{
    return general()(row,cols);
}

// column view (vector view)
template <typename FS>
const typename HeMatrix<FS>::ConstVectorView
HeMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col) const
{
    return general()(_,col);
}

template <typename FS>
typename HeMatrix<FS>::VectorView
HeMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col)
{
    return general()(_,col);
}

template <typename FS>
const typename HeMatrix<FS>::ConstVectorView
HeMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col) const
{
    return general()(rows,col);
}

template <typename FS>
typename HeMatrix<FS>::VectorView
HeMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col)
{
    return general()(rows,col);
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

// symmetric views
template <typename FS>
const typename HeMatrix<FS>::ConstSymmetricView
HeMatrix<FS>::symmetric() const
{
    return ConstSymmetricView(_engine, upLo());
}

template <typename FS>
typename HeMatrix<FS>::SymmetricView
HeMatrix<FS>::symmetric()
{
    return SymmetricView(_engine, upLo());
}

// triangular views
template <typename FS>
const typename HeMatrix<FS>::ConstTriangularView
HeMatrix<FS>::triangular() const
{
    if (upLo()==Upper) {
        return general().upper();
    }
    return general().lower();
}

template <typename FS>
typename HeMatrix<FS>::TriangularView
HeMatrix<FS>::triangular()
{
    if (upLo()==Upper) {
        return general().upper();
    }
    return general().lower();
}

// diag views
template <typename FS>
const typename HeMatrix<FS>::ConstVectorView
HeMatrix<FS>::diag(IndexType d) const
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(d>=0);
    } else {
        ASSERT(d<=0);
    }
#   endif

    return general().diag(d);
}

template <typename FS>
typename HeMatrix<FS>::VectorView
HeMatrix<FS>::diag(IndexType d)
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(d>=0);
    } else {
        ASSERT(d<=0);
    }
#   endif

    return general().diag(d);
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
HeMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::numCols() const
{
    return _engine.numCols();
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

template <typename FS>
StorageOrder
HeMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
template <typename RHS>
bool
HeMatrix<FS>::resize(const HeMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
HeMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(dim, dim, firstIndex, firstIndex, value);
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
