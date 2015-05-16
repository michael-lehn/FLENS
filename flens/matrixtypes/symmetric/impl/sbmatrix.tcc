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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SBMATRIX_TCC
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SBMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/symmetric/impl/sbmatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
SbMatrix<FS>::SbMatrix(StorageUpLo upLo)
    : upLo_(upLo)
{
}

template <typename FS>
SbMatrix<FS>::SbMatrix(IndexType dim, StorageUpLo upLo, IndexType numOffDiags,
                       IndexType firstIndex)
      : engine_(dim, dim,
                (upLo==Lower) ? numOffDiags : 0,
                (upLo==Upper) ? numOffDiags : 0,
                firstIndex),
        upLo_(upLo)
{
}

template <typename FS>
SbMatrix<FS>::SbMatrix(const Engine &engine, StorageUpLo upLo)
    : engine_(engine), upLo_(upLo)
{
}

template <typename FS>
SbMatrix<FS>::SbMatrix(const SbMatrix &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SbMatrix<FS>::SbMatrix(const SbMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SbMatrix<FS>::SbMatrix(SbMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS, class>
SbMatrix<FS>::SbMatrix(SbMatrix<RHS> &&rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SbMatrix<FS>::SbMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator=(const SbMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
SbMatrix<FS> &
SbMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
SbMatrix<FS> &
SbMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
SbMatrix<FS> &
SbMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator=(const ElementType &alpha)
{
    if (upLo_==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) = alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) = alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator+=(const ElementType &alpha)
{
    if (upLo_==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) += alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) += alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator-=(const ElementType &alpha)
{
    if (upLo_==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) -= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) -= alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator*=(const ElementType &alpha)
{
    if (upLo_==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) *= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) *= alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator/=(const ElementType &alpha)
{
    if (upLo_==Lower) {
        for (IndexType i=-numOffDiags(); i<=0; ++i)
            (*this).diag(i) /= alpha;
    } else {
        for (IndexType i=numOffDiags(); i>=0; --i)
            (*this).diag(i) /= alpha;
    }
    return *this;
}

template <typename FS>
const typename SbMatrix<FS>::ElementType &
SbMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
#   endif

    return engine_(row, col);
}

template <typename FS>
typename SbMatrix<FS>::ElementType &
SbMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo()==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
#   endif

    return engine_(row, col);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
const typename SbMatrix<FS>::ConstGeneralView
SbMatrix<FS>::general() const
{
    return ConstGeneralView(engine_);
}

template <typename FS>
typename SbMatrix<FS>::GeneralView
SbMatrix<FS>::general()
{
    return GeneralView(engine_);
}

// hermitian views
template <typename FS>
const typename SbMatrix<FS>::ConstHermitianView
SbMatrix<FS>::hermitian() const
{
    return ConstHermitianView(engine_, upLo_);
}

template <typename FS>
typename SbMatrix<FS>::HermitianView
SbMatrix<FS>::hermitian()
{
    return HermitianView(engine_, upLo_);
}

// triangular views
template <typename FS>
const typename SbMatrix<FS>::ConstTriangularView
SbMatrix<FS>::triangular() const
{
    return ConstTriangularView(engine_, upLo_);
}

template <typename FS>
typename SbMatrix<FS>::TriangularView
SbMatrix<FS>::triangular()
{
    return TriangularView(engine_, upLo_);
}

template <typename FS>
const typename SbMatrix<FS>::ConstVectorView
SbMatrix<FS>::diag(IndexType diag) const
{
    if ((upLo_==Upper)^(diag>0)) {
        ConstVectorView(engine_.viewDiag(-diag, engine_.firstIndex()));
    }

    return ConstVectorView(engine_.viewDiag(diag, engine_.firstIndex()));
}

template <typename FS>
typename SbMatrix<FS>::VectorView
SbMatrix<FS>::diag(IndexType diag)
{
    if ((upLo_==Upper)^(diag>0)) {
         return VectorView(engine_.viewDiag(-diag, engine_.firstIndex()));
    }

    return VectorView(engine_.viewDiag(diag, engine_.firstIndex()));
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::dim() const
{
    return engine_.dim();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::firstIndex() const
{
    return engine_.firstIndex();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::lastIndex() const
{
    return engine_.lastIndex();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numOffDiags() const
{
    return (upLo_==Upper) ? engine_.numSuperDiags()
                          : engine_.numSubDiags();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
SbMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
const typename SbMatrix<FS>::ElementType *
SbMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename SbMatrix<FS>::ElementType *
SbMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
template <typename RHS>
bool
SbMatrix<FS>::resize(const SbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
SbMatrix<FS>::resize(IndexType dim,
                     IndexType numOffDiags, IndexType firstIndex,
                     const ElementType &value)
{
    return engine_.resize(dim, dim,
                          (upLo_==Lower) ? numOffDiags : 0,
                          (upLo_==Upper) ? numOffDiags : 0,
                          firstIndex, value);
}

template <typename FS>
bool
SbMatrix<FS>::fill(const ElementType &value)
{
    return engine_.fill(value);
}

template <typename FS>
bool
SbMatrix<FS>::fillRandom()
{
    return engine_.fillRandom();
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename SbMatrix<FS>::Engine &
SbMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename SbMatrix<FS>::Engine &
SbMatrix<FS>::engine()
{
    return engine_;
}

template <typename FS>
StorageUpLo
SbMatrix<FS>::upLo() const
{
    return upLo_;
}

template <typename FS>
StorageUpLo &
SbMatrix<FS>::upLo()
{
    return upLo_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SBMATRIX_TCC
