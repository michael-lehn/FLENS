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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/triangular/impl/tbmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename FS>
TbMatrix<FS>::TbMatrix(IndexType dim, StorageUpLo upLo,
                       IndexType numOffDiags, Diag diag)
    : _engine(dim, dim, (upLo==Upper) ? 0 : numOffDiags,
                        (upLo==Upper) ? numOffDiags : 0),
      _upLo(upLo), _diag(diag)
{
}

template <typename FS>
TbMatrix<FS>::TbMatrix(const Engine &engine, StorageUpLo upLo, Diag diag)
    : _engine(engine), _upLo(upLo), _diag(diag)
{
}

template <typename FS>
TbMatrix<FS>::TbMatrix(const TbMatrix &rhs)
    : TriangularMatrix<TbMatrix<FS> >(),
      _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(const TbMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(TbMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const TbMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const ElementType &alpha)
{
    ASSERT(_diag!=NonUnit);

    if (_upLo==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) = alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) = alpha;
    }
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator+=(const ElementType &alpha)
{
    ASSERT(_diag!=NonUnit);

    if (_upLo==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) += alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) += alpha;
    }
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator-=(const ElementType &alpha)
{
    ASSERT(_diag!=NonUnit);

    if (_upLo==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) -= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) -= alpha;
    }
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator*=(const ElementType &alpha)
{
    ASSERT(_diag!=NonUnit);

    if (_upLo==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) *= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) *= alpha;
    }
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator/=(const ElementType &alpha)
{
    ASSERT(_diag!=NonUnit);

    if (_upLo==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) /= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) /= alpha;
    }
    return *this;
}

template <typename FS>
const typename TbMatrix<FS>::ElementType &
TbMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (_upLo==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((_diag==Unit) && (col==row)));
#   endif
    return _engine(row, col);
}

template <typename FS>
typename TbMatrix<FS>::ElementType &
TbMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (_upLo==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((_diag==Unit) && (col==row)));
#   endif
    return _engine(row, col);
}

// -- views --------------------------------------------------------------------
// general views
template <typename FS>
const typename TbMatrix<FS>::ConstGeneralView
TbMatrix<FS>::general() const
{
    return ConstGeneralView(_engine);
}

template <typename FS>
typename TbMatrix<FS>::GeneralView
TbMatrix<FS>::general()
{
    return GeneralView(_engine);
}

// hermitian views
template <typename FS>
const typename TbMatrix<FS>::ConstHermitianView
TbMatrix<FS>::hermitian() const
{
    ASSERT(_diag==NonUnit);
    return ConstHermitianView(_engine, _upLo);
}

template <typename FS>
typename TbMatrix<FS>::HermitianView
TbMatrix<FS>::hermitian()
{
    ASSERT(_diag==NonUnit);
    return HermitianView(_engine, _upLo);
}

// symmetric views
template <typename FS>
const typename TbMatrix<FS>::ConstSymmetricView
TbMatrix<FS>::symmetric() const
{
    ASSERT(_diag==NonUnit);
    return ConstSymmetricView(_engine, _upLo);
}

template <typename FS>
typename TbMatrix<FS>::SymmetricView
TbMatrix<FS>::symmetric()
{
    ASSERT(_diag==NonUnit);
    return SymmetricView(_engine, _upLo);
}

// triangular views
template <typename FS>
const typename TbMatrix<FS>::ConstView
TbMatrix<FS>::triangular() const
{
    return ConstView(_engine, _upLo);
}

template <typename FS>
typename TbMatrix<FS>::View
TbMatrix<FS>::triangular()
{
    return View(_engine, _upLo);
}


// diag views
template <typename FS>
const typename TbMatrix<FS>::ConstVectorView
TbMatrix<FS>::diag(IndexType diag) const
{
    return ConstVectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

template <typename FS>
typename TbMatrix<FS>::VectorView
TbMatrix<FS>::diag(IndexType diag)
{
    return VectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::dim() const
{
    ASSERT(_engine.numRows()==_engine.numCols());

    return _engine.numRows();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numCols() const
{
    return _engine.numCols();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::firstIndex() const
{
    ASSERT(_engine.firstRow()==_engine.firstCol());
    return _engine.firstRow();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numOffDiags() const
{
    return (_upLo==Upper)
             ? _engine.numSuperDiags()
             : _engine.numSubDiags();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numSubDiags() const
{
    return _engine.numSubDiags();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numSuperDiags() const
{
    return _engine.numSuperDiags();
}

template <typename FS>
const typename TbMatrix<FS>::ElementType *
TbMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename TbMatrix<FS>::ElementType *
TbMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
TbMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
template <typename RHS>
bool
TbMatrix<FS>::resize(const TbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
TbMatrix<FS>::resize(IndexType dim, IndexType numOffDiags,
                     IndexType firstIndex,
                     const ElementType &value)
{
    const IndexType numSubDiags = (_upLo == Upper) ? 0 : numOffDiags;
    const IndexType numSuperDiags = (_upLo == Upper) ? numOffDiags : 0;
    return _engine.resize(dim, dim, numSubDiags, numSuperDiags,
                          firstIndex, value);
}


template <typename FS>
bool
TbMatrix<FS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}

template <typename FS>
bool
TbMatrix<FS>::fillRandom()
{
    return _engine.fillRandom();
}

// -- implementation -----------------------------------------------------------
template <typename FS>
const typename TbMatrix<FS>::Engine &
TbMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename TbMatrix<FS>::Engine &
TbMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
TbMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
StorageUpLo &
TbMatrix<FS>::upLo()
{
    return _upLo;
}

template <typename FS>
Diag
TbMatrix<FS>::diag() const
{
    return _diag;
}

template <typename FS>
Diag &
TbMatrix<FS>::diag()
{
    return _diag;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC
