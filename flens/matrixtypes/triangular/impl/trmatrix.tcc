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

namespace flens {

template <typename FS>
TrMatrix<FS>::TrMatrix(const Engine &engine,
                       cxxblas::StorageUpLo upLo,
                       cxxblas::Diag diag)
    : _engine(engine), _upLo(upLo), _diag(diag)
{
}

template <typename FS>
TrMatrix<FS>::TrMatrix(const TrMatrix &rhs)
    : TriangularMatrix<TrMatrix<FS> >(),
      _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(const TrMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(TrMatrix<RHS> &rhs)
    : _engine(rhs.engine()), _upLo(rhs.upLo()), _diag(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TrMatrix<FS>::TrMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
const typename TrMatrix<FS>::ElementType &
TrMatrix<FS>::operator()(IndexType row, IndexType col) const
{
    ASSERT((_upLo==cxxblas::Upper) ? (col>=row) : (col<=row));
    ASSERT(!((_diag==cxxblas::Unit) && (col==row)));

    return _engine(row, col);
}

template <typename FS>
typename TrMatrix<FS>::ElementType &
TrMatrix<FS>::operator()(IndexType row, IndexType col)
{
    ASSERT((_upLo==cxxblas::Upper) ? (col>=row) : (col<=row));
    ASSERT(!((_diag==cxxblas::Unit) && (col==row)));

    return _engine(row, col);
}

// -- views ------------------------------------------------------------
// general views
template <typename FS>
typename TrMatrix<FS>::ConstGeneralView
TrMatrix<FS>::general() const
{
    return ConstGeneralView(_engine);
}

template <typename FS>
typename TrMatrix<FS>::GeneralView
TrMatrix<FS>::general()
{
    return GeneralView(_engine);
}

// hermitian views
template <typename FS>
typename TrMatrix<FS>::ConstHermitianView
TrMatrix<FS>::hermitian() const
{
    return ConstHermitianView(_engine, upLo());
}

template <typename FS>
typename TrMatrix<FS>::HermitianView
TrMatrix<FS>::hermitian()
{
    return HermitianView(_engine, upLo());
}

// symmetric views
template <typename FS>
typename TrMatrix<FS>::ConstSymmetricView
TrMatrix<FS>::symmetric() const
{
    ASSERT(diag()==cxxblas::NonUnit);
    return ConstSymmetricView(_engine, upLo());
}

template <typename FS>
typename TrMatrix<FS>::SymmetricView
TrMatrix<FS>::symmetric()
{
    ASSERT(diag()==cxxblas::NonUnit);
    return SymmetricView(_engine, upLo());
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::dim() const
{
    assert(_engine.numRows()==_engine.numCols());

    return _engine.numRows();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename FS>
typename TrMatrix<FS>::IndexType
TrMatrix<FS>::lastCol() const
{
    return _engine.lastCol();
}

// -- implementation -----------------------------------------------------------
template <typename FS>
const typename TrMatrix<FS>::Engine &
TrMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename TrMatrix<FS>::Engine &
TrMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
cxxblas::StorageUpLo
TrMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
cxxblas::StorageUpLo &
TrMatrix<FS>::upLo()
{
    return _upLo;
}

template <typename FS>
cxxblas::Diag
TrMatrix<FS>::diag() const
{
    return _diag;
}

template <typename FS>
cxxblas::Diag &
TrMatrix<FS>::diag()
{
    return _diag;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRMATRIX_TCC
