/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HECRSMATRIX_TCC
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HECRSMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/symmetric/impl/sycrsmatrix.h>

namespace flens {

// -- constructors -------------------------------------------------------------
template <typename CRS>
HeCRSMatrix<CRS>::HeCRSMatrix()
    : _upLo(Upper)
{
}

template <typename CRS>
template <typename RHS>
HeCRSMatrix<CRS>::HeCRSMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------
template <typename CRS>
template <typename RHS>
void
HeCRSMatrix<CRS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- methods ------------------------------------------------------------------
template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::dim() const
{
    ASSERT(numRows()==numCols());
    return _engine.numRows();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::numRows() const
{
    return _engine.numRows();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::numCols() const
{
    return _engine.numCols();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::indexBase() const
{
    return _engine.indexBase();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::firstRow() const
{
    return _engine.firstRow();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::lastRow() const
{
    return _engine.lastRow();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::firstCol() const
{
    return _engine.firstCol();
}

template <typename CRS>
typename HeCRSMatrix<CRS>::IndexType
HeCRSMatrix<CRS>::lastCol() const
{
    return _engine.lastCol();
}

template <typename CRS>
StorageUpLo
HeCRSMatrix<CRS>::upLo() const
{
    return _upLo;
}

template <typename CRS>
StorageUpLo &
HeCRSMatrix<CRS>::upLo()
{
    return _upLo;
}

// -- implementation -----------------------------------------------------------
template <typename CRS>
const typename HeCRSMatrix<CRS>::Engine &
HeCRSMatrix<CRS>::engine() const
{
    return _engine;
}

template <typename CRS>
typename HeCRSMatrix<CRS>::Engine &
HeCRSMatrix<CRS>::engine()
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HECRSMATRIX_TCC
