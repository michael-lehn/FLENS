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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCCSMATRIX_TCC
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCCSMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/symmetric/impl/syccsmatrix.h>

namespace flens {

// -- constructors -------------------------------------------------------------
template <typename CCS>
SyCCSMatrix<CCS>::SyCCSMatrix()
    : upLo_(Upper)
{
}

template <typename CCS>
template <typename RHS>
SyCCSMatrix<CCS>::SyCCSMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------
template <typename CCS>
template <typename RHS>
void
SyCCSMatrix<CCS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- methods ------------------------------------------------------------------
template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::dim() const
{
    ASSERT(numRows()==numCols());
    return engine_.numRows();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::numRows() const
{
    return engine_.numRows();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::indexBase() const
{
    return engine_.indexBase();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::numCols() const
{
    return engine_.numCols();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename CCS>
typename SyCCSMatrix<CCS>::IndexType
SyCCSMatrix<CCS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename CCS>
StorageUpLo
SyCCSMatrix<CCS>::upLo() const
{
    return upLo_;
}

template <typename CCS>
StorageUpLo &
SyCCSMatrix<CCS>::upLo()
{
    return upLo_;
}

// -- implementation -----------------------------------------------------------
template <typename CCS>
const typename SyCCSMatrix<CCS>::Engine &
SyCCSMatrix<CCS>::engine() const
{
    return engine_;
}

template <typename CCS>
typename SyCCSMatrix<CCS>::Engine &
SyCCSMatrix<CCS>::engine()
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCCSMATRIX_TCC
