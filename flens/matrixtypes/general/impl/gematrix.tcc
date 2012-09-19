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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/general/impl/ge/constelementclosure.tcc>
#include <flens/matrixtypes/general/impl/ge/elementclosure.tcc>
#include <flens/matrixtypes/general/impl/ge/initializer.tcc>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename FS>
GeMatrix<FS>::GeMatrix()
{
}

template <typename FS>
GeMatrix<FS>::GeMatrix(IndexType numRows, IndexType numCols)
    : _engine(numRows, numCols)
{
    ASSERT(numRows>=0);
    ASSERT(numCols>=0);
}

template <typename FS>
GeMatrix<FS>::GeMatrix(IndexType numRows, IndexType numCols,
                       IndexType firstRow, IndexType firstCol)
    : _engine(numRows, numCols, firstRow, firstCol)
{
    ASSERT(numRows>=0);
    ASSERT(numCols>=0);
}

template <typename FS>
GeMatrix<FS>::GeMatrix(const Range<IndexType> &rowRange,
                       const Range<IndexType> &colRange)
    : _engine(rowRange.numTicks(), colRange.numTicks(),
              rowRange.firstIndex(), colRange.firstIndex())
{
    ASSERT(rowRange.stride()==1);
    ASSERT(colRange.stride()==1);
}

template <typename FS>
GeMatrix<FS>::GeMatrix(const Engine &engine)
    : _engine(engine)
{
}

template <typename FS>
GeMatrix<FS>::GeMatrix(const GeMatrix &rhs)
    : GeneralMatrix<GeMatrix>(), _engine(rhs._engine)
{
}

template <typename FS>
template <typename RHS>
GeMatrix<FS>::GeMatrix(const GeMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GeMatrix<FS>::GeMatrix(GeMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GeMatrix<FS>::GeMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

template <typename FS>
template <typename VECTOR>
GeMatrix<FS>::GeMatrix(IndexType numRows, IndexType numCols, VECTOR &&rhs)
    : _engine(numRows, numCols, rhs.engine(), (FS::order==RowMajor) ? numCols
                                                                    : numRows)
{
}

template <typename FS>
template <typename VECTOR>
GeMatrix<FS>::GeMatrix(IndexType numRows, IndexType numCols,
                       VECTOR &&rhs,
                       IndexType leadingDimension)
    : _engine(numRows, numCols, rhs.engine(), leadingDimension)
{
}

// -- operators ----------------------------------------------------------------

template <typename FS>
typename GeMatrix<FS>::Initializer
GeMatrix<FS>::operator=(const ElementType &value)
{
    engine().fill(value);
    return Initializer(*this, firstRow(), firstCol());
}

template <typename FS>
GeMatrix<FS> &
GeMatrix<FS>::operator=(const GeMatrix<FS> &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
GeMatrix<FS> &
GeMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
GeMatrix<FS> &
GeMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
GeMatrix<FS> &
GeMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
GeMatrix<FS> &
GeMatrix<FS>::operator+=(const ElementType &alpha)
{
    const Underscore<IndexType> _;

    if (order()==ColMajor) {
        for (IndexType j=firstCol(); j<=lastCol(); ++j) {
            (*this)(_,j) += alpha;
        }
    } else {
        for (IndexType i=firstRow(); i<=lastRow(); ++i) {
            (*this)(i,_) += alpha;
        }
    }
    return *this;
}

template <typename FS>
GeMatrix<FS> &
GeMatrix<FS>::operator-=(const ElementType &alpha)
{
    const Underscore<IndexType> _;

    if (order()==ColMajor) {
        for (IndexType j=firstCol(); j<=lastCol(); ++j) {
            (*this)(_,j) -= alpha;
        }
    } else {
        for (IndexType i=firstRow(); i<=lastRow(); ++i) {
            (*this)(i,_) -= alpha;
        }
    }
    return *this;
}

template <typename FS>
GeMatrix<FS> &
GeMatrix<FS>::operator*=(const ElementType &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename FS>
GeMatrix<FS> &
GeMatrix<FS>::operator/=(const ElementType &alpha)
{
    blas::rscal(alpha, *this);
    return *this;
}

template <typename FS>
inline
const typename GeMatrix<FS>::ElementType &
GeMatrix<FS>::operator()(IndexType row, IndexType col) const
{
    return _engine(row, col);
}

template <typename FS>
inline
typename GeMatrix<FS>::ElementType &
GeMatrix<FS>::operator()(IndexType row, IndexType col)
{
    return _engine(row, col);
}

template <typename FS>
template <typename S>
const gematrix::ConstElementClosure<GeMatrix<FS>, typename Scalar<S>::Impl>
GeMatrix<FS>::operator()(const Scalar<S> &row, const Scalar<S> &col) const
{
    typedef typename Scalar<S>::Impl ScalarImpl;
    typedef gematrix::ConstElementClosure<GeMatrix, ScalarImpl>  CEC;
    return CEC(*this, row.impl(), col.impl());

}

template <typename FS>
const typename GeMatrix<FS>::ConstElementClosure
GeMatrix<FS>::operator()(const IndexVariable &row,
                         const IndexVariable &col) const
{
    return ConstElementClosure(*this, row, col);
}

template <typename FS>
typename GeMatrix<FS>::ElementClosure
GeMatrix<FS>::operator()(IndexVariable &row, IndexVariable &col)
{
    return ElementClosure(*this, row, col);
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::numCols() const
{
    return _engine.numCols();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::lastCol() const
{
    return _engine.lastCol();
}

template <typename FS>
Range<typename GeMatrix<FS>::IndexType>
GeMatrix<FS>::rows() const
{
    return Range<IndexType>(_engine.firstRow(),_engine.lastRow());
}

template <typename FS>
Range<typename GeMatrix<FS>::IndexType>
GeMatrix<FS>::cols() const
{
    return Range<IndexType>(_engine.firstCol(),_engine.lastCol());
}

template <typename FS>
const typename GeMatrix<FS>::ElementType *
GeMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename GeMatrix<FS>::ElementType *
GeMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
typename GeMatrix<FS>::IndexType
GeMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
GeMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
template <typename RHS>
bool
GeMatrix<FS>::resize(const GeMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
GeMatrix<FS>::resize(IndexType numRows, IndexType numCols,
                     IndexType firstRowIndex,
                     IndexType firstColIndex,
                     const ElementType &value)
{
    return _engine.resize(numRows, numCols,
                          firstRowIndex, firstColIndex,
                          value);
}

template <typename FS>
bool
GeMatrix<FS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}

template <typename FS>
void
GeMatrix<FS>::changeIndexBase(IndexType firstRowIndex, IndexType firstColIndex)
{
    _engine.changeIndexBase(firstRowIndex, firstColIndex);
}


// -- views --------------------------------------------------------------------
// vectorize matrix
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::vectorView() const
{
    return _engine.arrayView();
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::vectorView()
{
    return _engine.arrayView();
}

// vectorize matrix and select range
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::vectorView(IndexType from, IndexType to) const
{
    return _engine.arrayView().view(from,to);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::vectorView(IndexType from, IndexType to)
{
    return _engine.arrayView().view(from,to);
}

// diag views
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::diag(IndexType d) const
{
    return _engine.viewDiag(d);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::diag(IndexType d)
{
    return _engine.viewDiag(d);
}

// diag views
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::antiDiag(IndexType d) const
{
    return _engine.viewAntiDiag(d);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::antiDiag(IndexType d)
{
    return _engine.viewAntiDiag(d);
}

// triangular views
template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::upper() const
{
    return ConstTriangularView(engine(), Upper);
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::upper()
{
    return TriangularView(engine(), Upper);
}

template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::upperUnit() const
{
    return ConstTriangularView(engine(), Upper, Unit);
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::upperUnit()
{
    return TriangularView(engine(), Upper, Unit);
}

template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::strictUpper() const
{
    const Underscore<IndexType> _;
    const IndexType n = numCols();

    return operator()(_,_(2,n)).upper();
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::strictUpper()
{
    const Underscore<IndexType> _;
    const IndexType n = numCols();

    return operator()(_,_(2,n)).upper();
}

template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::lower() const
{
    return ConstTriangularView(engine(), Lower);
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::lower()
{
    return TriangularView(engine(), Lower);
}

template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::lowerUnit() const
{
    return ConstTriangularView(engine(), Lower, Unit);
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::lowerUnit()
{
    return TriangularView(engine(), Lower, Unit);
}

template <typename FS>
const typename GeMatrix<FS>::ConstTriangularView
GeMatrix<FS>::strictLower() const
{
    const Underscore<IndexType> _;
    const IndexType m = numRows();

    return operator()(_(2,m),_).lower();
}

template <typename FS>
typename GeMatrix<FS>::TriangularView
GeMatrix<FS>::strictLower()
{
    const Underscore<IndexType> _;
    const IndexType m = numRows();

    return operator()(_(2,m),_).lower();
}

// rectangular views
template <typename FS>
const typename GeMatrix<FS>::ConstView
GeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols) const
{
    ASSERT(rows.stride()==IndexType(1));
    ASSERT(cols.stride()==IndexType(1));
    return engine().view(rows.firstIndex(), cols.firstIndex(),
                         rows.lastIndex(), cols.lastIndex());
}

template <typename FS>
typename GeMatrix<FS>::View
GeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols)
{
    ASSERT(rows.stride()==IndexType(1));
    ASSERT(cols.stride()==IndexType(1));
    return engine().view(rows.firstIndex(), cols.firstIndex(),
                         rows.lastIndex(), cols.lastIndex());
}

template <typename FS>
template <typename RHS>
const typename GeMatrix<FS>::ConstView
GeMatrix<FS>::operator()(const GeMatrix<RHS> &A) const
{
    return engine().view(A.firstRow(), A.firstCol(),
                         A.lastRow(), A.lastCol());
}

template <typename FS>
template <typename RHS>
typename GeMatrix<FS>::View
GeMatrix<FS>::operator()(const GeMatrix<RHS> &A)
{
    return engine().view(A.firstRow(), A.firstCol(),
                         A.lastRow(), A.lastCol());
}

// rectangular views (all rows selected)
template <typename FS>
const typename GeMatrix<FS>::ConstView
GeMatrix<FS>::operator()(const Underscore<IndexType> &,
                         const Range<IndexType> &cols) const
{
    ASSERT(cols.stride()==IndexType(1));
    return engine().view(firstRow(), cols.firstIndex(),
                         lastRow(), cols.lastIndex());
}

template <typename FS>
typename GeMatrix<FS>::View
GeMatrix<FS>::operator()(const Underscore<IndexType> &,
                         const Range<IndexType> &cols)
{
    ASSERT(cols.stride()==IndexType(1));
    return engine().view(firstRow(), cols.firstIndex(),
                         lastRow(), cols.lastIndex());
}

// rectangular views (all columns selected)
template <typename FS>
const typename GeMatrix<FS>::ConstView
GeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &) const
{
    ASSERT(rows.stride()==IndexType(1));
    return engine().view(rows.firstIndex(), firstCol(),
                         rows.lastIndex(), lastCol());
}

template <typename FS>
typename GeMatrix<FS>::View
GeMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &)
{
    ASSERT(rows.stride()==IndexType(1));
    return engine().view(rows.firstIndex(), firstCol(),
                         rows.lastIndex(), lastCol());
}

// row view (vector view)
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &) const
{
    return engine().viewRow(row);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &)
{
    return engine().viewRow(row);
}

template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols) const
{
    return engine().viewRow(row, cols.firstIndex(), cols.lastIndex(), cols.stride());
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols)
{
    return engine().viewRow(row, cols.firstIndex(), cols.lastIndex(), cols.stride());
}

// column view (vector view)
template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::operator()(const Underscore<IndexType> &, IndexType col) const
{
    return engine().viewCol(col);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::operator()(const Underscore<IndexType> &, IndexType col)
{
    return engine().viewCol(col);
}

template <typename FS>
const typename GeMatrix<FS>::ConstVectorView
GeMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col) const
{
    return engine().viewCol(rows.firstIndex(), rows.lastIndex(), rows.stride(), col);
}

template <typename FS>
typename GeMatrix<FS>::VectorView
GeMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col)
{
    return engine().viewCol(rows.firstIndex(), rows.lastIndex(), rows.stride(), col);
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename GeMatrix<FS>::Engine &
GeMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename GeMatrix<FS>::Engine &
GeMatrix<FS>::engine()
{
    return _engine;
}

//-- GeMatrix specific functions -----------------------------------------------

//
//  imag
//
template <typename MZ>
ImagConstMatrixClosure<GeMatrix<MZ> >
imag(const GeMatrix<MZ> &Z)
{
    return Z;
}

template <typename MZ>
typename RestrictTo<IsGeMatrix<MZ>::value,
         ImagMatrixClosure<MZ> >::Type
imag(MZ &&Z)
{
    return Z;
}

//
//  real
//
template <typename MZ>
RealConstMatrixClosure<GeMatrix<MZ> >
real(const GeMatrix<MZ> &Z)
{
    return Z;
}

template <typename MZ>
typename RestrictTo<IsGeMatrix<MZ>::value,
         RealMatrixClosure<MZ> >::Type
real(MZ &&Z)
{
    return Z;
}

//
//  fillRandom
//

template <typename MA>
typename RestrictTo<IsGeMatrix<MA>::value,
         bool>::Type
fillRandom(MA &&A)
{
    fillRandom(A.engine());
}


} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_TCC
