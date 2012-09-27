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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/hermitian/impl/hpmatrix.h>
#include <flens/storage/storage.h>
#include <flens/typedefs.h>

namespace flens {

//-- Constructors --------------------------------------------------------------

template <typename PS>
HpMatrix<PS>::HpMatrix(IndexType dim, StorageUpLo upLo)
    : _engine(dim), _upLo(upLo)
{
    ASSERT(dim>=0);
}


template <typename PS>
HpMatrix<PS>::HpMatrix(const Engine &engine, StorageUpLo upLo)
    : _engine(engine), _upLo(upLo)
{
}

template <typename PS>
HpMatrix<PS>::HpMatrix(const HpMatrix &rhs)
    : HermitianMatrix<HpMatrix<PS> >(),
      _engine(rhs.engine()), _upLo(upLo)
{
}

template <typename PS>
template <typename RHS>
HpMatrix<PS>::HpMatrix(const HpMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(upLo)
{
}

template <typename PS>
template <typename RHS>
HpMatrix<PS>::HpMatrix(HpMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(upLo)
{
}

template <typename PS>
template <typename RHS>
HpMatrix<PS>::HpMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- Operators -----------------------------------------------------------------

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    VectorView x = ArrayView(_engine.numNonZeros(), _engine.data());

    x = alpha;
    return *this;
}

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator=(const HpMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename PS>
template <typename RHS>
HpMatrix<PS> &
HpMatrix<PS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename PS>
template <typename RHS>
HpMatrix<PS> &
HpMatrix<PS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename PS>
template <typename RHS>
HpMatrix<PS> &
HpMatrix<PS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename PS>
const typename HpMatrix<PS>::ElementType &
HpMatrix<PS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
#   endif
    return _engine(upLo(), row, col);
}

template <typename PS>
typename HpMatrix<PS>::ElementType &
HpMatrix<PS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
#   endif
    return _engine(upLo(), row, col);
}

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator+=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    VectorView x = ArrayView(_engine.numNonZeros(), _engine.data());

    x += alpha;
    return *this;
}

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator-=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    VectorView x = ArrayView(_engine.numNonZeros(), _engine.data());

    x -= alpha;
    return *this;
}

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator*=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    VectorView x = ArrayView(_engine.numNonZeros(), _engine.data());

    x *= alpha;
    return *this;
}

template <typename PS>
HpMatrix<PS> &
HpMatrix<PS>::operator/=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    VectorView x = ArrayView(_engine.numNonZeros(), _engine.data());

    x /= alpha;
    return *this;
}

//-- Methods -------------------------------------------------------------------

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::numRows() const
{
    return _engine.dim();
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::numCols() const
{
    return _engine.dim();
}


template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::dim() const
{
    return _engine.dim();
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::firstRow() const
{
    return _engine.indexBase();
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::lastRow() const
{
    return firstRow()+numRows()-1;
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::firstCol() const
{
    return _engine.indexBase();
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::lastCol() const
{
    return firstCol()+numCols()-1;
}

template <typename PS>
typename HpMatrix<PS>::IndexType
HpMatrix<PS>::indexBase() const
{
    return _engine.indexBase();
}

template <typename PS>
StorageUpLo
HpMatrix<PS>::upLo() const
{
    return _upLo;
}

template <typename PS>
StorageUpLo &
HpMatrix<PS>::upLo()
{
    return _upLo;
}

template <typename PS>
const typename HpMatrix<PS>::ElementType *
HpMatrix<PS>::data() const
{
    return _engine.data();
}

template <typename PS>
typename HpMatrix<PS>::ElementType *
HpMatrix<PS>::data()
{
    return _engine.data();
}

template <typename PS>
StorageOrder
HpMatrix<PS>::order() const
{
    return _engine.order;
}

template <typename PS>
bool
HpMatrix<PS>::fill(const ElementType &value)
{
    ASSERT(cxxblas::imag(value)==0);

    return _engine.fill(value);
}

template <typename PS>
template <typename RHS>
bool
HpMatrix<PS>::resize(const HpMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.dim(), rhs.indexBase(), value);
}

template <typename PS>
bool
HpMatrix<PS>::resize(IndexType dim, IndexType indexBase,
                     const ElementType &value)
{
    return _engine.resize(dim, indexBase, value);
}

// -- views --------------------------------------------------------------------

// general views
template <typename PS>
const typename HpMatrix<PS>::ConstView
HpMatrix<PS>::hermitian() const
{
    return ConstView(_engine);
}

template <typename PS>
typename HpMatrix<PS>::View
HpMatrix<PS>::hermitian()
{
    return View(_engine);
}

// symmetric view
template <typename PS>
const typename HpMatrix<PS>::ConstSymmetricView
HpMatrix<PS>::symmetric() const
{
    return ConstSymmetricView(_engine);
}

template <typename PS>
typename HpMatrix<PS>::SymmetricView
HpMatrix<PS>::symmetric()
{
    return SymmetricView(_engine);
}

// triangular view
template <typename PS>
const typename HpMatrix<PS>::ConstTriangularView
HpMatrix<PS>::triangular() const
{
    return ConstTriangularView(_engine);
}

template <typename PS>
typename HpMatrix<PS>::TriangularView
HpMatrix<PS>::triangular()
{
    return TriangularView(_engine);
}

// -- implementation -----------------------------------------------------------

template <typename PS>
const typename HpMatrix<PS>::Engine &
HpMatrix<PS>::engine() const
{
    return _engine;
}

template <typename PS>
typename HpMatrix<PS>::Engine &
HpMatrix<PS>::engine()
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC
