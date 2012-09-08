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

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
GbMatrix<FS>::GbMatrix()
{
}

template <typename FS>
GbMatrix<FS>::GbMatrix(IndexType n, IndexType m, IndexType numSubDiags, IndexType numSuperDiags, IndexType firstIndex)
    : _engine(n, m, numSubDiags, numSuperDiags, firstIndex)
{
    ASSERT(n>=0);
    ASSERT(m>=0);
    ASSERT(numSubDiags>=0);
    ASSERT(numSuperDiags>=0);
}

template <typename FS>
GbMatrix<FS>::GbMatrix(const Engine &engine)
    : _engine(engine)
{
}

template <typename FS>
GbMatrix<FS>::GbMatrix(const GbMatrix &rhs)
    : GeneralMatrix<GbMatrix<FS> >(),
      _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GbMatrix<FS>::GbMatrix(const GbMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
GbMatrix<FS>::GbMatrix(GbMatrix<RHS> &rhs)
    : _engine(rhs.engine())
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
        (*this).viewDiag(i) = alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator+=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        (*this).viewDiag(i) += alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator-=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        (*this).viewDiag(i) -= alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator*=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        (*this).viewDiag(i) *= alpha;
    return *this;
}

template <typename FS>
GbMatrix<FS> &
GbMatrix<FS>::operator/=(const ElementType &alpha)
{
    for (IndexType i = -numSubDiags(); i <= numSuperDiags(); ++i)
        (*this).viewDiag(i) /= alpha;
    return *this;
}

template <typename FS>
const typename GbMatrix<FS>::ElementType &
GbMatrix<FS>::operator()(IndexType row, IndexType col) const
{
    return _engine(row, col);
}

template <typename FS>
typename GbMatrix<FS>::ElementType &
GbMatrix<FS>::operator()(IndexType row, IndexType col)
{
    return _engine(row, col);
}

// -- views --------------------------------------------------------------------

// general views


// triangular view
template <typename FS>
const typename GbMatrix<FS>::ConstTriangularView
GbMatrix<FS>::upper(Diag diag) const
{
    ASSERT(_engine.numCols() == _engine.numRows());
    return ConstTriangularView(_engine.viewDiags(0, _engine.numSuperDiags()), Upper, diag);
}

template <typename FS>
typename GbMatrix<FS>::TriangularView
GbMatrix<FS>::upper(Diag diag)
{
    ASSERT(_engine.numCols() == _engine.numRows());
    return TriangularView(_engine.viewDiags(0, _engine.numSuperDiags()), Upper, diag);
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::strictUpper() const
{
    ASSERT(_engine.numSuperDiags()>0);
    return ConstView(_engine.viewDiags(1, _engine.numSuperDiags()));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::strictUpper() 
{
    ASSERT(_engine.numSuperDiags()>0);
    return View(_engine.viewDiags(1, _engine.numSuperDiags()));
}

template <typename FS>
const typename GbMatrix<FS>::ConstTriangularView
GbMatrix<FS>::lower(Diag diag) const
{
    ASSERT(_engine.numCols() == _engine.numRows());
    return ConstTriangularView(_engine.viewDiags(-_engine.numSubDiags(), 0), Lower, diag);
}

template <typename FS>
typename GbMatrix<FS>::TriangularView
GbMatrix<FS>::lower(Diag diag)
{
    ASSERT(_engine.numCols() == _engine.numRows());
    return TriangularView(_engine.viewDiags(-_engine.numSubDiags(), 0), Lower, diag);
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::strictLower() const
{
    ASSERT(_engine.numSubDiags()>0);
    return ConstView(_engine.viewDiags(-_engine.numSubDiags(),-1));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::strictLower() 
{
    ASSERT(_engine.numSubDiags()>0);
    return View(_engine.viewDiags(-_engine.numSubDiags(),-1));
}

template <typename FS>
const typename GbMatrix<FS>::ConstVectorView
GbMatrix<FS>::viewDiag(IndexType diag) const
{
    return ConstVectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

template <typename FS>
typename GbMatrix<FS>::VectorView
GbMatrix<FS>::viewDiag(IndexType diag)
{
    return VectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

template <typename FS>
const typename GbMatrix<FS>::ConstView
GbMatrix<FS>::viewDiags(IndexType fromDiag, IndexType toDiag) const
{
    return ConstView(_engine.viewDiags(fromDiag, toDiag));
}

template <typename FS>
typename GbMatrix<FS>::View
GbMatrix<FS>::viewDiags(IndexType fromDiag, IndexType toDiag) 
{
    return View(_engine.viewDiags(fromDiag, toDiag));
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numCols() const
{
    return _engine.numCols();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastCol() const
{
    return _engine.lastCol();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numSubDiags() const
{
    return _engine.numSubDiags();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::numSuperDiags() const
{
    return _engine.numSuperDiags();
}

template <typename FS>
typename GbMatrix<FS>::IndexType
GbMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
GbMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
const typename GbMatrix<FS>::ElementType *
GbMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename GbMatrix<FS>::ElementType *
GbMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
template <typename RHS>
bool
GbMatrix<FS>::resize(const GbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
GbMatrix<FS>::resize(IndexType n, IndexType m, IndexType kl, IndexType ku, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(n, m, kl, ku, firstIndex, value);
}

template <typename FS>
bool
GbMatrix<FS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}

template <typename FS>
bool
GbMatrix<FS>::fillRandom()
{
    return _engine.fillRandom();
}
// -- implementation -----------------------------------------------------------

template <typename FS>
const typename GbMatrix<FS>::Engine &
GbMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename GbMatrix<FS>::Engine &
GbMatrix<FS>::engine()
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GBMATRIX_TCC
