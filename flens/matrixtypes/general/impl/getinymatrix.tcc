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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/general/impl/getinymatrix.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename TFS>
GeTinyMatrix<TFS>::GeTinyMatrix()
{
}

template <typename TFS>
GeTinyMatrix<TFS>::GeTinyMatrix(const Engine &engine)
    : engine_(engine)
{
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS>::GeTinyMatrix(const GeTinyMatrix<RHS> &rhs)
    : engine_(rhs.engine)
{
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS>::GeTinyMatrix(GeTinyMatrix<RHS> &rhs)
    : engine_(rhs.engine)
{
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS>::GeTinyMatrix(const Matrix<RHS> &rhs)
{
    assign(NoTrans, rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename TFS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator=(const GeTinyMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename TFS>
template <typename RHS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename TFS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator+=(const ElementType &alpha)
{
    for (IndexType i=firstRow(); i<=lastRow(); ++i) {
        for (IndexType j=firstCol(); j<=lastCol(); ++j) {
            (*this)(i,j) += alpha;
        }
    }
    return *this;
}

template <typename TFS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator-=(const ElementType &alpha)
{
    for (IndexType i=firstRow(); i<=lastRow(); ++i) {
        for (IndexType j=firstCol(); j<=lastCol(); ++j) {
            (*this)(i,j) -= alpha;
        }
    }
    return *this;
}

template <typename TFS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator*=(const ElementType &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename TFS>
GeTinyMatrix<TFS> &
GeTinyMatrix<TFS>::operator/=(const ElementType &alpha)
{
    blas::rscal(alpha, *this);
    return *this;
}

template <typename TFS>
const typename GeTinyMatrix<TFS>::ElementType &
GeTinyMatrix<TFS>::operator()(IndexType row, IndexType col) const
{
    return engine_(row, col);
}

template <typename TFS>
typename GeTinyMatrix<TFS>::ElementType &
GeTinyMatrix<TFS>::operator()(IndexType row, IndexType col)
{
    return engine_(row, col);
}

// -- methods ------------------------------------------------------------------
template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::numRows() const
{
    return Engine::numRows;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::numCols() const
{
    return Engine::numCols;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::firstRow() const
{
    return Engine::firstRow;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::lastRow() const
{
    return Engine::lastRow;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::firstCol() const
{
    return Engine::firstCol;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::lastCol() const
{
    return Engine::lastCol;
}

template <typename TFS>
const typename GeTinyMatrix<TFS>::ElementType *
GeTinyMatrix<TFS>::data() const
{
    return engine_.data();
}

template <typename TFS>
typename GeTinyMatrix<TFS>::ElementType *
GeTinyMatrix<TFS>::data()
{
    return engine_.data();
}

template <typename TFS>
typename GeTinyMatrix<TFS>::IndexType
GeTinyMatrix<TFS>::leadingDimension() const
{
    return Engine::leadingDimension;
}

template <typename TFS>
void
GeTinyMatrix<TFS>::fill(const ElementType &value)
{
    engine_.fill(value);
}

// -- implementation -----------------------------------------------------------
template <typename TFS>
const typename GeTinyMatrix<TFS>::Engine &
GeTinyMatrix<TFS>::engine() const
{
    return engine_;
}

template <typename TFS>
typename GeTinyMatrix<TFS>::Engine &
GeTinyMatrix<TFS>::engine()
{
    return engine_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_TCC
