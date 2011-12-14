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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_TCC
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_TCC 1

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/symmetric/impl/symatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
SyMatrix<FS>::SyMatrix(const Engine &engine, StorageUpLo upLo)
    : _engine(engine), _upLo(upLo)
{
    ASSERT(_engine.numRows()==_engine.numCols());
}

template <typename FS>
SyMatrix<FS>::SyMatrix(const SyMatrix &rhs)
    : SymmetricMatrix<SyMatrix<FS> >(),
      _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(const SyMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(SyMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SyMatrix<FS>::SyMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
void
SyMatrix<FS>::operator=(const ElementType &value)
{
    engine().fill(_upLo, value);
}

template <typename FS>
SyMatrix<FS> &
SyMatrix<FS>::operator=(const SyMatrix &rhs)
{
    blas::copy(NoTrans, rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
SyMatrix<FS> &
SyMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    blas::copy(NoTrans, rhs, *this);
    return *this;
}

template <typename FS>
const typename SyMatrix<FS>::ElementType &
SyMatrix<FS>::operator()(IndexType row, IndexType col) const
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
typename SyMatrix<FS>::ElementType &
SyMatrix<FS>::operator()(IndexType row, IndexType col)
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
const typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols) const
{
    return general()(rows, cols);
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols)
{
    return general()(rows, cols);
}

// rectangular views (all rows selected)
template <typename FS>
const typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols) const
{
    return general()(_,cols);
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols)
{
    return general()(_,cols);
}

// rectangular views (all columns selected)
template <typename FS>
const typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_) const
{
    return general()(rows,_);
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_)
{
    return general()(rows,_);
}

// row view (vector view)
template <typename FS>
const typename SyMatrix<FS>::ConstVectorView
SyMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_) const
{
    return general()(row,_);
}

template <typename FS>
typename SyMatrix<FS>::VectorView
SyMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_)
{
    return general()(row,_);
}

template <typename FS>
const typename SyMatrix<FS>::ConstVectorView
SyMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols) const
{
    return general()(row,cols);
}

template <typename FS>
typename SyMatrix<FS>::VectorView
SyMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols)
{
    return general()(row,cols);
}

// column view (vector view)
template <typename FS>
const typename SyMatrix<FS>::ConstVectorView
SyMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col) const
{
    return general()(_,col);
}

template <typename FS>
typename SyMatrix<FS>::VectorView
SyMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col)
{
    return general()(_,col);
}

template <typename FS>
const typename SyMatrix<FS>::ConstVectorView
SyMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col) const
{
    return general()(rows,col);
}

template <typename FS>
typename SyMatrix<FS>::VectorView
SyMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col)
{
    return general()(rows,col);
}

// -- views --------------------------------------------------------------------
// general views
template <typename FS>
const typename SyMatrix<FS>::ConstGeneralView
SyMatrix<FS>::general() const
{
    return ConstGeneralView(_engine);
}

template <typename FS>
typename SyMatrix<FS>::GeneralView
SyMatrix<FS>::general()
{
    return GeneralView(_engine);
}

// hermitian views
template <typename FS>
const typename SyMatrix<FS>::ConstHermitianView
SyMatrix<FS>::hermitian() const
{
    return ConstHermitianView(_engine, upLo());
}

template <typename FS>
typename SyMatrix<FS>::HermitianView
SyMatrix<FS>::hermitian()
{
    return HermitianView(_engine, upLo());
}

// triangular views
template <typename FS>
const typename SyMatrix<FS>::ConstTriangularView
SyMatrix<FS>::triangular() const
{
    if (upLo()==Upper) {
        return general().upper();
    }
    return general().lower();
}

template <typename FS>
typename SyMatrix<FS>::TriangularView
SyMatrix<FS>::triangular()
{
    if (upLo()==Upper) {
        return general().upper();
    }
    return general().lower();
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::dim() const
{
    ASSERT(_engine.numRows()==_engine.numCols());

    return _engine.numRows();
}

template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::lastCol() const
{
    return _engine.lastCol();
}

template <typename FS>
const typename SyMatrix<FS>::ElementType *
SyMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename SyMatrix<FS>::ElementType *
SyMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
typename SyMatrix<FS>::IndexType
SyMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
SyMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
template <typename RHS>
bool
SyMatrix<FS>::resize(const SyMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
SyMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(dim, dim, firstIndex, firstIndex, value);
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename SyMatrix<FS>::Engine &
SyMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename SyMatrix<FS>::Engine &
SyMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
SyMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
StorageUpLo &
SyMatrix<FS>::upLo()
{
    return _upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_TCC
