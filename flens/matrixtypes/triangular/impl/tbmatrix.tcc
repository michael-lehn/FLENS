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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/triangular/impl/tbmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename FS>
TbMatrix<FS>::TbMatrix(IndexType dim, StorageUpLo upLo,
                       IndexType numOffDiags, Diag diag)
    : engine_(dim, dim, (upLo==Upper) ? 0 : numOffDiags,
                        (upLo==Upper) ? numOffDiags : 0),
      upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TbMatrix<FS>::TbMatrix(const Engine &engine, StorageUpLo upLo, Diag diag)
    : engine_(engine), upLo_(upLo), diag_(diag)
{
}

template <typename FS>
TbMatrix<FS>::TbMatrix(const TbMatrix &rhs)
    : TriangularMatrix<TbMatrix<FS> >(),
      engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(const TbMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(TbMatrix<RHS> &rhs)
    : engine_(rhs.engine()), upLo_(rhs.upLo()), diag_(rhs.diag())
{
}

template <typename FS>
template <typename RHS>
TbMatrix<FS>::TbMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const TbMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
TbMatrix<FS> &
TbMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
TbMatrix<FS> &
TbMatrix<FS>::operator=(const ElementType &alpha)
{
    ASSERT(diag_!=NonUnit);

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
TbMatrix<FS> &
TbMatrix<FS>::operator+=(const ElementType &alpha)
{
    ASSERT(diag_!=NonUnit);

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
TbMatrix<FS> &
TbMatrix<FS>::operator-=(const ElementType &alpha)
{
    ASSERT(diag_!=NonUnit);

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
TbMatrix<FS> &
TbMatrix<FS>::operator*=(const ElementType &alpha)
{
    ASSERT(diag_!=NonUnit);

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
TbMatrix<FS> &
TbMatrix<FS>::operator/=(const ElementType &alpha)
{
    ASSERT(diag_!=NonUnit);

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
const typename TbMatrix<FS>::ElementType &
TbMatrix<FS>::operator()(IndexType row, IndexType col) const
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((diag_==Unit) && (col==row)));
#   endif
    return engine_(row, col);
}

template <typename FS>
typename TbMatrix<FS>::ElementType &
TbMatrix<FS>::operator()(IndexType row, IndexType col)
{
#   ifndef NDEBUG
    if (upLo_==Upper) {
        ASSERT(col>=row);
    } else {
        ASSERT(col<=row);
    }
    ASSERT(!((diag_==Unit) && (col==row)));
#   endif
    return engine_(row, col);
}

// -- views --------------------------------------------------------------------
// general views
template <typename FS>
const typename TbMatrix<FS>::ConstGeneralView
TbMatrix<FS>::general() const
{
    return ConstGeneralView(engine_);
}

template <typename FS>
typename TbMatrix<FS>::GeneralView
TbMatrix<FS>::general()
{
    return GeneralView(engine_);
}

// hermitian views
template <typename FS>
const typename TbMatrix<FS>::ConstHermitianView
TbMatrix<FS>::hermitian() const
{
    ASSERT(diag_==NonUnit);
    return ConstHermitianView(engine_, upLo_);
}

template <typename FS>
typename TbMatrix<FS>::HermitianView
TbMatrix<FS>::hermitian()
{
    ASSERT(diag_==NonUnit);
    return HermitianView(engine_, upLo_);
}

// symmetric views
template <typename FS>
const typename TbMatrix<FS>::ConstSymmetricView
TbMatrix<FS>::symmetric() const
{
    ASSERT(diag_==NonUnit);
    return ConstSymmetricView(engine_, upLo_);
}

template <typename FS>
typename TbMatrix<FS>::SymmetricView
TbMatrix<FS>::symmetric()
{
    ASSERT(diag_==NonUnit);
    return SymmetricView(engine_, upLo_);
}

// triangular views
template <typename FS>
const typename TbMatrix<FS>::ConstView
TbMatrix<FS>::triangular() const
{
    return ConstView(engine_, upLo_);
}

template <typename FS>
typename TbMatrix<FS>::View
TbMatrix<FS>::triangular()
{
    return View(engine_, upLo_);
}


// diag views
template <typename FS>
const typename TbMatrix<FS>::ConstVectorView
TbMatrix<FS>::diag(IndexType diag) const
{
    return ConstVectorView(engine_.viewDiag(diag, engine_.firstIndex()));
}

template <typename FS>
typename TbMatrix<FS>::VectorView
TbMatrix<FS>::diag(IndexType diag)
{
    return VectorView(engine_.viewDiag(diag, engine_.firstIndex()));
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::dim() const
{
    ASSERT(engine_.numRows()==engine_.numCols());

    return engine_.numRows();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::firstIndex() const
{
    ASSERT(engine_.firstRow()==engine_.firstCol());
    return engine_.firstRow();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::lastIndex() const
{
    return engine_.lastIndex();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numOffDiags() const
{
    return (upLo_==Upper)
             ? engine_.numSuperDiags()
             : engine_.numSubDiags();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numSubDiags() const
{
    return engine_.numSubDiags();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::numSuperDiags() const
{
    return engine_.numSuperDiags();
}

template <typename FS>
const typename TbMatrix<FS>::ElementType *
TbMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename TbMatrix<FS>::ElementType *
TbMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
typename TbMatrix<FS>::IndexType
TbMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
TbMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
template <typename RHS>
bool
TbMatrix<FS>::resize(const TbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
TbMatrix<FS>::resize(IndexType dim, IndexType numOffDiags,
                     IndexType firstIndex,
                     const ElementType &value)
{
    const IndexType numSubDiags = (upLo_ == Upper) ? 0 : numOffDiags;
    const IndexType numSuperDiags = (upLo_ == Upper) ? numOffDiags : 0;
    return engine_.resize(dim, dim, numSubDiags, numSuperDiags,
                          firstIndex, value);
}


template <typename FS>
bool
TbMatrix<FS>::fill(const ElementType &value)
{
    return engine_.fill(value);
}

template <typename FS>
bool
TbMatrix<FS>::fillRandom()
{
    return engine_.fillRandom();
}

// -- implementation -----------------------------------------------------------
template <typename FS>
const typename TbMatrix<FS>::Engine &
TbMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename TbMatrix<FS>::Engine &
TbMatrix<FS>::engine()
{
    return engine_;
}

template <typename FS>
StorageUpLo
TbMatrix<FS>::upLo() const
{
    return upLo_;
}

template <typename FS>
StorageUpLo &
TbMatrix<FS>::upLo()
{
    return upLo_;
}

template <typename FS>
Diag
TbMatrix<FS>::diag() const
{
    return diag_;
}

template <typename FS>
Diag &
TbMatrix<FS>::diag()
{
    return diag_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_TCC
