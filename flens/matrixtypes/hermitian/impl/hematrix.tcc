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
#include <flens/matrixtypes/hermitian/impl/hematrix.h>
#include <flens/typedefs.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename FS>
HeMatrix<FS>::HeMatrix(StorageUpLo upLo)
    : upLo_(upLo)
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(IndexType dim, StorageUpLo upLo)
    : engine_(dim, dim), upLo_(upLo)
{
    ASSERT(dim>=0);
}

template <typename FS>
HeMatrix<FS>::HeMatrix(IndexType dim, IndexType firstIndex, StorageUpLo upLo)
    : engine_(dim, dim, firstIndex, firstIndex), upLo_(upLo)
{
    ASSERT(dim>=0);
}

template <typename FS>
HeMatrix<FS>::HeMatrix(IndexType dim, IndexType firstRow, IndexType firstCol,
                       StorageUpLo upLo)
    : engine_(dim, dim, firstRow, firstCol), upLo_(upLo)
{
    ASSERT(dim>=0);
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const Engine &engine, StorageUpLo upLo)
    : engine_(engine), upLo_(upLo)
{
    ASSERT(engine_.numRows()==engine_.numCols());
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const HeMatrix &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const HeMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(HeMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS, class>
HeMatrix<FS>::HeMatrix(HeMatrix<RHS> &&rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
void
HeMatrix<FS>::operator=(const ElementType &value)
{
    engine().fill(upLo_, value);
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
template <typename RHS>
HeMatrix<FS> &
HeMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
HeMatrix<FS> &
HeMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
const typename HeMatrix<FS>::ElementType &
HeMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
#   endif
    return engine_(row, col);
}

template <typename FS>
typename HeMatrix<FS>::ElementType &
HeMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
#   endif
    return engine_(row, col);
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

// -- methods ------------------------------------------------------------------

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::dim() const
{
    ASSERT(engine_.numRows()==engine_.numCols());

    return engine_.numRows();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename FS>
const typename HeMatrix<FS>::ElementType *
HeMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename HeMatrix<FS>::ElementType *
HeMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
typename HeMatrix<FS>::IndexType
HeMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
HeMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
bool
HeMatrix<FS>::fill(const ElementType &value)
{
    ASSERT(cxxblas::imag(value)==0);
    return engine_.fill(upLo_, value);
}

template <typename FS>
template <typename RHS>
bool
HeMatrix<FS>::resize(const HeMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
HeMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                     const ElementType &value)
{
    return engine_.resize(dim, dim, firstIndex, firstIndex, value);
}

template <typename FS>
bool
HeMatrix<FS>::resize(IndexType dim, StorageUpLo upLo, IndexType firstIndex,
                     const ElementType &value)
{
    upLo_ = upLo;
    return engine_.resize(dim, dim, firstIndex, firstIndex, value);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
const typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::general() const
{
    return ConstGeneralView(engine_);
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::general()
{
    return GeneralView(engine_);
}

// symmetric views
template <typename FS>
const typename HeMatrix<FS>::ConstSymmetricView
HeMatrix<FS>::symmetric() const
{
    return ConstSymmetricView(engine_, upLo());
}

template <typename FS>
typename HeMatrix<FS>::SymmetricView
HeMatrix<FS>::symmetric()
{
    return SymmetricView(engine_, upLo());
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


// -- implementation -----------------------------------------------------------

template <typename FS>
const typename HeMatrix<FS>::Engine &
HeMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename HeMatrix<FS>::Engine &
HeMatrix<FS>::engine()
{
    return engine_;
}

template <typename FS>
StorageUpLo
HeMatrix<FS>::upLo() const
{
    return upLo_;
}

template <typename FS>
StorageUpLo &
HeMatrix<FS>::upLo()
{
    return upLo_;
}

//-- HeMatrix specific functions -----------------------------------------------

//
//  fillRandom
//

template <typename MA>
typename RestrictTo<IsHeMatrix<MA>::value,
         bool>::Type
fillRandom(MA &&A)
{
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::ElementType   ElementType;
    typedef typename MatrixA::IndexType     IndexType;
    typedef typename MatrixA::VectorView    VectorView;

    fillRandom(A.upLo(), A.engine());
    VectorView d = A.diag(0);
    for (IndexType i=d.firstIndex();i<=d.lastIndex();++i) {
        d(i) = ElementType(cxxblas::real(d(i)));
    }
    return true;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HEMATRIX_TCC
