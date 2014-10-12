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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRMATRIX_TCC 1

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/symmetric/impl/symatrix.h>
#include <flens/matrixtypes/triangular/impl/trmatrix.h>
#include <flens/matrixtypes/triangular/impl/tr/constelementclosure.tcc>
#include <flens/matrixtypes/triangular/impl/tr/elementclosure.tcc>
#include <flens/typedefs.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename FS>
TrMatrix<FS>::TrMatrix(IndexType dim, StorageUpLo upLo, Diag diag)
    : engine_(dim, dim), upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TrMatrix<FS>::TrMatrix(IndexType numRows, IndexType numCols,
                       StorageUpLo upLo, Diag diag)
    : engine_(numRows, numCols), upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TrMatrix<FS>::TrMatrix(IndexType numRows, IndexType numCols,
                       IndexType firstRow, IndexType firstCol,
                       StorageUpLo upLo, Diag diag)
    : engine_(numRows, numCols, firstRow, firstCol), upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TrMatrix<FS>::TrMatrix(const Engine &engine, StorageUpLo upLo, Diag diag)
    : engine_(engine), upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TrMatrix<FS>::TrMatrix(const TrMatrix &rhs)
    : TriangularMatrix<TrMatrix<FS> >(),
      engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(const TrMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(TrMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------
template <typename FS>
void
TrMatrix<FS>::operator=(const ElementType &value)
{
    engine().fill(upLo_, value);
}

template <typename FS>
TrMatrix<FS> &
TrMatrix<FS>::operator=(const TrMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
TrMatrix<FS> &
TrMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TrMatrix<FS> &
TrMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    ASSERT(this->diag()==NonUnit);
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TrMatrix<FS> &
TrMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    ASSERT(this->diag()==NonUnit);
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
const typename TrMatrix<FS>::ElementType &
TrMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
    ASSERT(!((diag_==Unit) && (col==row)));
#   endif
    return engine_(row, col);
}

template <typename FS>
typename TrMatrix<FS>::ElementType &
TrMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
    ASSERT(!((diag_==Unit) && (col==row)));
#   endif
    return engine_(row, col);
}

template <typename FS>
template <typename S>
const trmatrix::ConstElementClosure<TrMatrix<FS>, typename Scalar<S>::Impl>
TrMatrix<FS>::operator()(const Scalar<S> &row, const Scalar<S> &col) const
{
    typedef typename Scalar<S>::Impl ScalarImpl;
    typedef trmatrix::ConstElementClosure<TrMatrix, ScalarImpl>  CEC;
    return CEC(*this, row.impl(), col.impl());

}

template <typename FS>
const typename TrMatrix<FS>::ConstElementClosure
TrMatrix<FS>::operator()(const IndexVariable &row,
                         const IndexVariable &col) const
{
    return ConstElementClosure(*this, row, col);
}

template <typename FS>
typename TrMatrix<FS>::ElementClosure
TrMatrix<FS>::operator()(IndexVariable &row, IndexVariable &col)
{
    return ElementClosure(*this, row, col);
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::dim() const
{
    ASSERT(engine_.numRows()==engine_.numCols());

    return engine_.numRows();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename FS>
StorageUpLo
TrMatrix<FS>::upLo() const
{
    return upLo_;
}

template <typename FS>
StorageUpLo &
TrMatrix<FS>::upLo()
{
    return upLo_;
}

template <typename FS>
Diag
TrMatrix<FS>::diag() const
{
    return diag_;
}

template <typename FS>
Diag &
TrMatrix<FS>::diag()
{
    return diag_;
}

template <typename FS>
const typename TrMatrix<FS>::ElementType *
TrMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename TrMatrix<FS>::ElementType *
TrMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
TrMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
bool
TrMatrix<FS>::fill(const ElementType &value)
{
    return engine_.fill(upLo_, value);
}

template <typename FS>
template <typename RHS>
bool
TrMatrix<FS>::resize(const TrMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
TrMatrix<FS>::resize(IndexType numRows, IndexType numCols,
                     IndexType firstRowIndex, IndexType firstColIndex,
                     const ElementType &value)
{
    return engine_.resize(numRows, numCols,
                          firstRowIndex, firstColIndex,
                          value);
}

// -- views ------------------------------------------------------------
// diag views
template <typename FS>
const typename TrMatrix<FS>::ConstVectorView
TrMatrix<FS>::diag(IndexType d) const
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
typename TrMatrix<FS>::VectorView
TrMatrix<FS>::diag(IndexType d)
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

// general views
template <typename FS>
const typename TrMatrix<FS>::ConstGeneralView
TrMatrix<FS>::general() const
{
    return ConstGeneralView(engine_);
}

template <typename FS>
typename TrMatrix<FS>::GeneralView
TrMatrix<FS>::general()
{
    return GeneralView(engine_);
}

// hermitian views
template <typename FS>
const typename TrMatrix<FS>::ConstHermitianView
TrMatrix<FS>::hermitian() const
{
    ASSERT(numRows()==numCols());
    return ConstHermitianView(engine_, upLo());
}

template <typename FS>
typename TrMatrix<FS>::HermitianView
TrMatrix<FS>::hermitian()
{
    ASSERT(numRows()==numCols());
    return HermitianView(engine_, upLo());
}

// symmetric views
template <typename FS>
const typename TrMatrix<FS>::ConstSymmetricView
TrMatrix<FS>::symmetric() const
{
    ASSERT(numRows()==numCols());
    return ConstSymmetricView(engine_, upLo());
}

template <typename FS>
typename TrMatrix<FS>::SymmetricView
TrMatrix<FS>::symmetric()
{
    ASSERT(numRows()==numCols());
    return SymmetricView(engine_, upLo());
}

// rectangular views
template <typename FS>
const typename TrMatrix<FS>::ConstGeneralView
TrMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols) const
{
    return general()(rows, cols);
}

template <typename FS>
typename TrMatrix<FS>::GeneralView
TrMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Range<IndexType> &cols)
{
    return general()(rows, cols);
}

// rectangular views (all rows selected)
template <typename FS>
const typename TrMatrix<FS>::ConstGeneralView
TrMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols) const
{
    return general()(_,cols);
}

template <typename FS>
typename TrMatrix<FS>::GeneralView
TrMatrix<FS>::operator()(const Underscore<IndexType> &_,
                         const Range<IndexType> &cols)
{
    return general()(_,cols);
}

// rectangular views (all columns selected)
template <typename FS>
const typename TrMatrix<FS>::ConstGeneralView
TrMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_) const
{
    return general()(rows,_);
}

template <typename FS>
typename TrMatrix<FS>::GeneralView
TrMatrix<FS>::operator()(const Range<IndexType> &rows,
                         const Underscore<IndexType> &_)
{
    return general()(rows,_);
}

// row view (vector view)
template <typename FS>
const typename TrMatrix<FS>::ConstVectorView
TrMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_) const
{
    return general()(row,_);
}

template <typename FS>
typename TrMatrix<FS>::VectorView
TrMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &_)
{
    return general()(row,_);
}

template <typename FS>
const typename TrMatrix<FS>::ConstVectorView
TrMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols) const
{
    return general()(row,cols);
}

template <typename FS>
typename TrMatrix<FS>::VectorView
TrMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols)
{
    return general()(row,cols);
}

// column view (vector view)
template <typename FS>
const typename TrMatrix<FS>::ConstVectorView
TrMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col) const
{
    return general()(_,col);
}

template <typename FS>
typename TrMatrix<FS>::VectorView
TrMatrix<FS>::operator()(const Underscore<IndexType> &_, IndexType col)
{
    return general()(_,col);
}

template <typename FS>
const typename TrMatrix<FS>::ConstVectorView
TrMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col) const
{
    return general()(rows,col);
}

template <typename FS>
typename TrMatrix<FS>::VectorView
TrMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col)
{
    return general()(rows,col);
}

// -- implementation -----------------------------------------------------------
template <typename FS>
const typename TrMatrix<FS>::Engine &
TrMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename TrMatrix<FS>::Engine &
TrMatrix<FS>::engine()
{
    return engine_;
}

//-- TrMatrix specific functions -----------------------------------------------

//
//  fillRandom
//

template <typename MA>
typename RestrictTo<IsTrMatrix<MA>::value,
         bool>::Type
fillRandom(MA &&A)
{
    fillRandom(A.upLo(), A.engine());
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRMATRIX_TCC
