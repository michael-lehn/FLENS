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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GECRSMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GECRSMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/general/impl/gecrsmatrix.h>

namespace flens {

// -- constructor --------------------------------------------------------------
template <typename CRS>
GeCRSMatrix<CRS>::GeCRSMatrix()
{
}

template <typename CRS>
template <typename RHS>
GeCRSMatrix<CRS>::GeCRSMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- operators -----------------------------------------------------------------
template <typename CRS>
template <typename RHS>
void
GeCRSMatrix<CRS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- methods -------------------------------------------------------------------
template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::numRows() const
{
    return engine_.numRows();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::numCols() const
{
    return engine_.numCols();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::indexBase() const
{
    return engine_.indexBase();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename CRS>
typename GeCRSMatrix<CRS>::IndexType
GeCRSMatrix<CRS>::lastCol() const
{
    return engine_.lastCol();
}

//-- implementation ------------------------------------------------------------
template <typename CRS>
const typename GeCRSMatrix<CRS>::Engine &
GeCRSMatrix<CRS>::engine() const
{
    return engine_;
}

template <typename CRS>
typename GeCRSMatrix<CRS>::Engine &
GeCRSMatrix<CRS>::engine()
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GECRSMATRIX_TCC

