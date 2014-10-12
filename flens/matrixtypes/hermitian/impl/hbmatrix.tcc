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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HBMATRIX_TCC
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HBMATRIX_TCC 1

#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/hermitian/impl/hbmatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
HbMatrix<FS>::HbMatrix()
{
}

template <typename FS>
HbMatrix<FS>::HbMatrix(IndexType dim, StorageUpLo upLo, IndexType numOffDiags,
                       IndexType firstIndex)
    : engine_(dim, dim, (upLo == Lower) ? numOffDiags : 0,
              (upLo == Upper) ? numOffDiags : 0, firstIndex),
      upLo_(upLo)
{
    ASSERT(dim>=0);
}

template <typename FS>
HbMatrix<FS>::HbMatrix(const Engine &engine, StorageUpLo upLo)
    : engine_(engine), upLo_(upLo)
{
}

template <typename FS>
HbMatrix<FS>::HbMatrix(const HbMatrix &rhs)
    : HermitianMatrix<HbMatrix<FS> >(),
      engine_(rhs.engine()),
      upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HbMatrix<FS>::HbMatrix(const HbMatrix<RHS> &rhs)
    : engine_(rhs.engine()),
      upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HbMatrix<FS>::HbMatrix(HbMatrix<RHS> &rhs)
    : engine_(rhs.engine()),
      upLo_(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HbMatrix<FS>::HbMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator=(const HbMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
HbMatrix<FS> &
HbMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
HbMatrix<FS> &
HbMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
HbMatrix<FS> &
HbMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    if (upLo_==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).diag(i) = alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).diag(i) = alpha;
    }
    return *this;
}

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator+=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    if (upLo_==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).diag(i) += alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).diag(i) += alpha;
    }
    return *this;
}

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator-=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    if (upLo_==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).diag(i) -= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).diag(i) -= alpha;
    }
    return *this;
}

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator*=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    if (upLo_==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).diag(i) *= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).diag(i) *= alpha;
    }
    return *this;
}

template <typename FS>
HbMatrix<FS> &
HbMatrix<FS>::operator/=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);

    if (upLo_==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).diag(i) /= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).diag(i) /= alpha;
    }
    return *this;
}

template <typename FS>
const typename HbMatrix<FS>::ElementType &
HbMatrix<FS>::operator()(IndexType row, IndexType col) const
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
typename HbMatrix<FS>::ElementType &
HbMatrix<FS>::operator()(IndexType row, IndexType col)
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
typename HbMatrix<FS>::ConstGeneralView
HbMatrix<FS>::general() const
{
    return ConstGeneralViewView(engine_, upLo_);
}

template <typename FS>
typename HbMatrix<FS>::GeneralView
HbMatrix<FS>::general()
{
    return GeneralView(engine_, upLo_);
}

// symmetric view
template <typename FS>
typename HbMatrix<FS>::ConstSymmetricView
HbMatrix<FS>::symmetric() const
{
    return ConstSymmetricView(engine_, upLo_);
}

template <typename FS>
typename HbMatrix<FS>::SymmetricView
HbMatrix<FS>::symmetric()
{
    return SymmetricView(engine_, upLo_);
}

// triangular view
template <typename FS>
typename HbMatrix<FS>::ConstTriangularView
HbMatrix<FS>::triangular() const
{
    return ConstTriangularView(engine_, upLo_);
}

template <typename FS>
typename HbMatrix<FS>::TriangularView
HbMatrix<FS>::triangular()
{
    return TriangularView(engine_, upLo_);
}

template <typename FS>
const typename HbMatrix<FS>::ConstVectorView
HbMatrix<FS>::diag(IndexType diag) const
{
    return ConstVectorView(engine_.viewDiag(diag));
}

template <typename FS>
typename HbMatrix<FS>::VectorView
HbMatrix<FS>::diag(IndexType diag)
{
    return VectorView(engine_.viewDiag(diag));
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::dim() const
{
    return engine_.dim();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numCols() const
{
    return engine_.numCols();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numRows() const
{
    return engine_.numRows();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::firstIndex() const
{
    return engine_.firstIndex();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::lastIndex() const
{
    return engine_.lastIndex();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::leadingDimension() const
{
    return engine_.leadingDimension();
}

template <typename FS>
StorageOrder
HbMatrix<FS>::order() const
{
    return engine_.order;
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numOffDiags() const
{
    return (upLo_==Upper) ? engine_.numSuperDiags()
                          : engine_.numSubDiags();
}

template <typename FS>
const typename HbMatrix<FS>::ElementType *
HbMatrix<FS>::data() const
{
    return engine_.data();
}

template <typename FS>
typename HbMatrix<FS>::ElementType *
HbMatrix<FS>::data()
{
    return engine_.data();
}

template <typename FS>
template <typename RHS>
bool
HbMatrix<FS>::resize(const HbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return engine_.resize(rhs.engine(), value);
}

template <typename FS>
bool
HbMatrix<FS>::resize(IndexType dim, IndexType numOffDiags,
                      IndexType firstIndex, const ElementType &value)
{
    return engine_.resize(dim, dim,
                          (upLo_==Lower) ? numOffDiags : 0,
                          (upLo_==Upper) ? numOffDiags : 0,
                          firstIndex, value);
}

template <typename FS>
bool
HbMatrix<FS>::fill(const ElementType &value)
{
    ASSERT(cxxblas::imag(value)==0);

    return engine_.fill(value);
}

template <typename FS>
bool
HbMatrix<FS>::fillRandom()
{
    bool val = engine_.fillRandom();
    VectorView d = (*this).diag(0);
    for (IndexType i=d.firstIndex();i<=d.lastIndex();++i) {
        d(i) = ElementType(cxxblas::real(d(i)));
    }
    return val;
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename HbMatrix<FS>::Engine &
HbMatrix<FS>::engine() const
{
    return engine_;
}

template <typename FS>
typename HbMatrix<FS>::Engine &
HbMatrix<FS>::engine()
{
    return engine_;
}

template <typename FS>
StorageUpLo
HbMatrix<FS>::upLo() const
{
    return upLo_;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HBMATRIX_TCC
