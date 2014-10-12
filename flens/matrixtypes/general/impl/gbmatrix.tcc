/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GBMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GBMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>
#include <flens/matrixtypes/general/impl/gematrix.h>

namespace flens {

template <typename FS>
GbMatrix<FS>::GbMatrix()
{
}

template <typename FS>
GbMatrix<FS>::GbMatrix(IndexType n, IndexType m,
                       IndexType numSubDiags,
                       IndexType numSuperDiags,
                       IndexType firstIndex)
    : engine_(n, m, numSubDiags, numSuperDiags, firstIndex)
{
    ASSERT(n>=0);
    ASSERT(m>=0);
    ASSERT(numSubDiags>=0);
    ASSERT(numSuperDiags>=0);
}

template <typename FS>
GbMatrix<FS>::GbMatrix(const Engine &engine)
    : engine_(engine)
{
}

template <typename FS>
GbMatrix<FS>::GbMatrix(const GbMatrix &rhs)
    : GeneralMatrix<GbMatrix<FS> >(),
      engine_(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GbMatrix<FS>::GbMatrix(const GbMatrix<RHS> &rhs)
    : engine_(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GbMatrix<FS>::GbMatrix(GbMatrix<RHS> &rhs)
    : engine_(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GbMatrix<FS>::GbMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------
template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator=(const GbMatrix<FS> &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
GbMatrix<FS> &
GbMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
GbMatrix<FS> &
GbMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
GbMatrix<FS> &
GbMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        diag(i) = alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator+=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        diag(i) += alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator-=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        diag(i) -= alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator*=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        diag(i) *= alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator/=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        diag(i) /= alpha;
    return *this;
}

template <typename FS>
const typename GbMatrix<FS>::ElementType &
GbMatrix<FS>::operator()(IndexType row, IndexType col) const
{
    return engine_(row, col);
}

template <typename FS>
typename GbMatrix<FS>::ElementType &
GbMatrix<FS>::operator()(IndexType row, IndexType col)
{
    return engine_(row, col);
}

// -- views --------------------------------------------------------------------

// general views


// triangular view
template <typename FS>
const typename GbMatrix<FS>::ConstTriangularView
GbMatrix<FS>::upper(Diag diag) const
{
    ASSERT(engine_.numCols() == engine_.numRows());
    return ConstTriangularView(engine_.viewDiags(0, engine_.numSuperDiags()),
                               Upper, diag);
}

template <typename FS>
typename GbMatrix<FS>::TriangularView
GbMatrix<FS>::upper(Diag diag)
{
    ASSERT(engine_.numCols() == engine_.numRows());
    return TriangularView(engine_.viewDiags(0, engine_.numSuperDiags()),
                          Upper, diag);
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::strictUpper() const
{
    ASSERT(engine_.numSuperDiags()>0);
    return ConstView(engine_.viewDiags(1, engine_.numSuperDiags()));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::strictUpper()
{
    ASSERT(engine_.numSuperDiags()>0);
    return View(engine_.viewDiags(1, engine_.numSuperDiags()));
}

template <typename FS>
const typename GbMatrix<FS>::ConstTriangularView
GbMatrix<FS>::lower(Diag diag) const
{
    ASSERT(engine_.numCols() == engine_.numRows());
    return ConstTriangularView(engine_.viewDiags(-engine_.numSubDiags(), 0),
                               Lower, diag);
}

template <typename FS>
typename GbMatrix<FS>::TriangularView
GbMatrix<FS>::lower(Diag diag)
{
    ASSERT(engine_.numCols() == engine_.numRows());
    return TriangularView(engine_.viewDiags(-engine_.numSubDiags(), 0),
                          Lower, diag);
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::strictLower() const
{
    ASSERT(engine_.numSubDiags()>0);
    return ConstView(engine_.viewDiags(-engine_.numSubDiags(),-1));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::strictLower()
{
    ASSERT(engine_.numSubDiags()>0);
    return View(engine_.viewDiags(-engine_.numSubDiags(),-1));
}

template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::diag(IndexType d) const
{
    return ConstVectorView(engine_.viewDiag(d, engine_.firstIndex()));
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::diag(IndexType d)
{
    return VectorView(engine_.viewDiag(d, engine_.firstIndex()));
}

// row view (vector view)
template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &) const
{
    return engine().viewRow(row);
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::operator()(IndexType row, const Underscore<IndexType> &)
{
    return engine().viewRow(row);
}

template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols) const
{
    const IndexType firstIndex = cols.firstIndex();
    const IndexType lastIndex  = cols.lastIndex();
    const IndexType stride     = cols.stride();

    return engine().viewRow(row, firstIndex, lastIndex, stride);
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::operator()(IndexType row, const Range<IndexType> &cols)
{
    const IndexType firstIndex = cols.firstIndex();
    const IndexType lastIndex  = cols.lastIndex();
    const IndexType stride     = cols.stride();

    return engine().viewRow(row, firstIndex, lastIndex, stride);
}

// column view (vector view)
template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::operator()(const Underscore<IndexType> &, IndexType col) const
{
    return engine().viewCol(col);
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::operator()(const Underscore<IndexType> &, IndexType col)
{
    return engine().viewCol(col);
}

template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col) const
{
    const IndexType firstIndex = rows.firstIndex();
    const IndexType lastIndex  = rows.lastIndex();
    const IndexType stride     = rows.stride();

    return engine().viewCol(firstIndex, lastIndex, stride, col);
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::operator()(const Range<IndexType> &rows, IndexType col)
{
    const IndexType firstIndex = rows.firstIndex();
    const IndexType lastIndex  = rows.lastIndex();
    const IndexType stride     = rows.stride();

    return engine().viewCol(firstIndex, lastIndex, stride, col);
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::diags(IndexType fromDiag, IndexType toDiag) const
{
    return ConstView(engine_.viewDiags(fromDiag, toDiag));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::diags(IndexType fromDiag, IndexType toDiag)
{
    return View(engine_.viewDiags(fromDiag, toDiag));
}

template <typename FS>
const typename GbMatrix<FS>::ConstGeView
GbMatrix<FS>::viewStorageGeMatrix() const
{
    return engine_.viewFullStorage();
}

template <typename FS>
typename GbMatrix<FS>::GeView
GbMatrix<FS>::viewStorageGeMatrix()
{
    return engine_.viewFullStorage();
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstIndex() const
{
    return engine_.firstIndex();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastIndex() const
{
    return engine_.lastIndex();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numSubDiags() const
{
    return engine_.numSubDiags();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numSuperDiags() const
{
    return engine_.numSuperDiags();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
GbMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
const typename GbMatrix<FS>::ElementType *
GbMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename GbMatrix<FS>::ElementType *
GbMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
template <typename RHS>
bool
GbMatrix<FS>::resize(const GbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
GbMatrix<FS>::resize(IndexType n, IndexType m, IndexType kl, IndexType ku,
                     IndexType firstIndex,
                     const ElementType &value)
{
    return engine_.resize(n, m, kl, ku, firstIndex, value);
}

template <typename FS>
bool
GbMatrix<FS>::fill(const ElementType &value)
{
    return engine_.fill(value);
}

template <typename FS>
bool
GbMatrix<FS>::fillRandom()
{
    return engine_.fillRandom();
}
// -- implementation -----------------------------------------------------------

template <typename FS>
const typename GbMatrix<FS>::Engine &
GbMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename GbMatrix<FS>::Engine &
GbMatrix<FS>::engine()
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GBMATRIX_TCC
