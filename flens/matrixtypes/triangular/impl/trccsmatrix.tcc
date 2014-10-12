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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCCSMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCCSMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/triangular/impl/trccsmatrix.h>

namespace flens {

// -- constructors -------------------------------------------------------------
template <typename CCS>
TrCCSMatrix<CCS>::TrCCSMatrix()
    : upLo_(Upper)
{
}

template <typename CCS>
template <typename RHS>
TrCCSMatrix<CCS>::TrCCSMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------
template <typename CCS>
template <typename RHS>
void
TrCCSMatrix<CCS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- methods ------------------------------------------------------------------
template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::dim() const
{
    ASSERT(numRows()==numCols());
    return engine_.numRows();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::numRows() const
{
    return engine_.numRows();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::indexBase() const
{
    return engine_.indexBase();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::numCols() const
{
    return engine_.numCols();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename CCS>
typename TrCCSMatrix<CCS>::IndexType
TrCCSMatrix<CCS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename CCS>
StorageUpLo
TrCCSMatrix<CCS>::upLo() const
{
    return upLo_;
}

template <typename CCS>
StorageUpLo &
TrCCSMatrix<CCS>::upLo()
{
    return upLo_;
}

// -- implementation -----------------------------------------------------------
template <typename CCS>
const typename TrCCSMatrix<CCS>::Engine &
TrCCSMatrix<CCS>::engine() const
{
    return engine_;
}

template <typename CCS>
typename TrCCSMatrix<CCS>::Engine &
TrCCSMatrix<CCS>::engine()
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCCSMATRIX_TCC
