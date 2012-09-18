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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/triangular/impl/tpmatrix.h>
#include <flens/typedefs.h>

namespace flens {

//
//-- Constructors --------------------------------------------------------------
//

template <typename PS>
TpMatrix<PS>::TpMatrix(IndexType dim, StorageUpLo upLo, Diag diag)
    : _engine(dim), _upLo(upLo), _diag(diag)
{
}


template <typename PS>
TpMatrix<PS>::TpMatrix(const Engine &engine, StorageUpLo upLo, Diag diag)
    : _engine(engine), _upLo(upLo), _diag(diag)
{
}

template <typename PS>
TpMatrix<PS>::TpMatrix(const TpMatrix &rhs)
    : TriangularMatrix<TpMatrix<PS> >(),
      _engine(rhs.engine()), _upLo(upLo), _diag(rhs.diag())
{
}

template <typename PS>
template <typename RHS>
TpMatrix<PS>::TpMatrix(const TpMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(upLo), _diag(rhs.diag())
{
}

template <typename PS>
template <typename RHS>
TpMatrix<PS>::TpMatrix(TpMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(upLo), _diag(rhs.diag())
{
}

template <typename PS>
template <typename RHS>
TpMatrix<PS>::TpMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- Operators -----------------------------------------------------------------

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator=(const ElementType &alpha)
{
    ASSERT(_diag==NonUnit);

    VectorView x = ArrayView(engine.numNonZeros(), engine.data());

    x = alpha;
    return *this;
}

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator=(const TpMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename PS>
template <typename RHS>
TpMatrix<PS> &
TpMatrix<PS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename PS>
const typename TpMatrix<PS>::ElementType &
TpMatrix<PS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((_diag==Unit) && (col==row)));
#   endif
    return _engine(upLo(), row, col);
}

template <typename PS>
typename TpMatrix<PS>::ElementType &
TpMatrix<PS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((_diag==Unit) && (col==row)));
#   endif
    return _engine(upLo(), row, col);
}

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator+=(const ElementType &alpha)
{
    ASSERT(_diag==NonUnit);

    VectorView x = ArrayView(engine.numNonZeros(), engine.data());

    x += alpha;
    return *this;
}

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator-=(const ElementType &alpha)
{
    ASSERT(_diag==NonUnit);

    VectorView x = ArrayView(engine.numNonZeros(), engine.data());

    x -= alpha;
    return *this;
}

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator*=(const ElementType &alpha)
{
    ASSERT(_diag==NonUnit);

    VectorView x = ArrayView(engine.numNonZeros(), engine.data());

    x *= alpha;
    return *this;
}

template <typename PS>
TpMatrix<PS> &
TpMatrix<PS>::operator/=(const ElementType &alpha)
{
    ASSERT(_diag==NonUnit);

    VectorView x = ArrayView(engine.numNonZeros(), engine.data());

    x /= alpha;
    return *this;
}

//-- Methods -------------------------------------------------------------------
template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::numRows() const
{
    return _engine.dim();
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::numCols() const
{
    return _engine.dim();
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::dim() const
{
    return _engine.dim();
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::firstRow() const
{
    return _engine.indexBase();
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::lastRow() const
{
    return firstRow()+numRows()-1;
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::firstCol() const
{
    return _engine.indexBase();
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::lastCol() const
{
    return firstCol()+numCols()-1;
}

template <typename PS>
typename TpMatrix<PS>::IndexType
TpMatrix<PS>::indexBase() const
{
    return _engine.indexBase();
}


template <typename PS>
StorageUpLo
TpMatrix<PS>::upLo() const
{
    return _upLo;
}

template <typename PS>
StorageUpLo &
TpMatrix<PS>::upLo()
{
    return _upLo;
}

template <typename PS>
Diag
TpMatrix<PS>::diag() const
{
    return _diag;
}

template <typename PS>
Diag &
TpMatrix<PS>::diag()
{
    return _diag;
}

template <typename PS>
const typename TpMatrix<PS>::ElementType *
TpMatrix<PS>::data() const
{
    return _engine.data();
}

template <typename PS>
typename TpMatrix<PS>::ElementType *
TpMatrix<PS>::data()
{
    return _engine.data();
}

template <typename PS>
StorageOrder
TpMatrix<PS>::order() const
{
    return _engine.order;
}

template <typename PS>
bool
TpMatrix<PS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}

template <typename PS>
template <typename RHS>
bool
TpMatrix<PS>::resize(const TpMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return resize(rhs.dim(), rhs.indexBase(), value);
}

template <typename PS>
bool
TpMatrix<PS>::resize(IndexType dim,
                     IndexType indexBase,
                     const ElementType &value)
{
    return _engine.resize(dim, indexBase,
                          value);
}

//-- Views ---------------------------------------------------------------------

// hermitian views
template <typename PS>
const typename TpMatrix<PS>::ConstHermitianView
TpMatrix<PS>::hermitian() const
{
    return ConstHermitianView(_engine, upLo());
}

template <typename PS>
typename TpMatrix<PS>::HermitianView
TpMatrix<PS>::hermitian()
{
    return HermitianView(_engine, upLo());
}

// symmetric views
template <typename PS>
const typename TpMatrix<PS>::ConstSymmetricView
TpMatrix<PS>::symmetric() const
{
    return ConstSymmetricView(_engine, upLo());
}

template <typename PS>
typename TpMatrix<PS>::SymmetricView
TpMatrix<PS>::symmetric()
{
    return SymmetricView(_engine, upLo());
}

// triangular views
template <typename PS>
const typename TpMatrix<PS>::ConstView
TpMatrix<PS>::triangular() const
{
    return ConstView(_engine, upLo(), diag());
}

template <typename PS>
typename TpMatrix<PS>::View
TpMatrix<PS>::triangular()
{
    return View(_engine, upLo(), diag());
}


//-- Implementation ------------------------------------------------------------

template <typename PS>
const typename TpMatrix<PS>::Engine &
TpMatrix<PS>::engine() const
{
    return _engine;
}

template <typename PS>
typename TpMatrix<PS>::Engine &
TpMatrix<PS>::engine()
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_TCC
