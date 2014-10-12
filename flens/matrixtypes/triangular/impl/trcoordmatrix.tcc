/*
 *   Copyright (c) 2012-2013, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCOORDMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCOORDMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/matrixtypes/triangular/impl/trcoordmatrix.h>

namespace flens {

template <typename CS>
TrCoordMatrix<CS>::TrCoordMatrix(IndexType   dim,
                                 StorageUpLo upLo,
                                 IndexType   densityEstimate,
                                 IndexType   indexBase)
    : engine_(dim, dim, densityEstimate, indexBase), upLo_(upLo)
{
}

// -- operators ----------------------------------------------------------------
template <typename CS>
typename TrCoordMatrix<CS>::ElementProxy
TrCoordMatrix<CS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col-firstCol()>=row-firstRow());
    } else {
        ASSERT(col-firstCol()<=row-firstRow());
    }
#   endif

    return engine_(row, col);
}

// -- methods ------------------------------------------------------------------
template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::dim() const
{
    ASSERT(numRows()==numCols());
    return engine_.numRows();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::numRows() const
{
    return engine_.numRows();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::numCols() const
{
    return engine_.numCols();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::indexBase() const
{
    return engine_.indexBase();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::firstRow() const
{
    return engine_.firstRow();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::lastRow() const
{
    return engine_.lastRow();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::firstCol() const
{
    return engine_.firstCol();
}

template <typename CS>
typename TrCoordMatrix<CS>::IndexType
TrCoordMatrix<CS>::lastCol() const
{
    return engine_.lastCol();
}

template <typename CS>
StorageUpLo
TrCoordMatrix<CS>::upLo() const
{
    return upLo_;
}

// -- implementation -----------------------------------------------------------
template <typename CS>
const typename TrCoordMatrix<CS>::Engine &
TrCoordMatrix<CS>::engine() const
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TRCOORDMATRIX_TCC
