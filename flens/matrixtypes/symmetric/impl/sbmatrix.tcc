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

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/symmetric/impl/symatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
SbMatrix<FS>::SbMatrix()
{
  
}
  
template <typename FS>
SbMatrix<FS>::SbMatrix(IndexType dim, StorageUpLo upLo, IndexType numOffDiags,
                       IndexType firstIndex)
      : _engine(dim, dim, (upLo == Lower) ? numOffDiags : 0, (upLo == Upper) ? numOffDiags : 0,
                firstIndex),
      _upLo(upLo)
{
  
}
  
template <typename FS>
SbMatrix<FS>::SbMatrix(const Engine &engine, StorageUpLo upLo)
    : _engine(engine), _upLo(upLo)
{
}

template <typename FS>
SbMatrix<FS>::SbMatrix(const SbMatrix &rhs)
    : SymmetricMatrix<SbMatrix<FS> >(),
      _engine(rhs.engine()),
      _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SbMatrix<FS>::SbMatrix(const SbMatrix<RHS> &rhs)
    : _engine(rhs.engine()),
      _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
SbMatrix<FS>::SbMatrix(SbMatrix<RHS> &rhs)
    : _engine(rhs.engine()),
      _upLo(rhs.upLo())
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
    if (_upLo==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).viewDiag(i) = alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).viewDiag(i) = alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator+=(const ElementType &alpha)
{
    if (_upLo==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).viewDiag(i) += alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).viewDiag(i) += alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator-=(const ElementType &alpha)
{
    if (_upLo==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).viewDiag(i) -= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).viewDiag(i) -= alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator*=(const ElementType &alpha)
{
    if (_upLo==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).viewDiag(i) *= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).viewDiag(i) *= alpha;
    }
    return *this;
}

template <typename FS>
SbMatrix<FS> &
SbMatrix<FS>::operator/=(const ElementType &alpha)
{
    if (_upLo==Lower) {
        for (IndexType i = -numOffDiags(); i <= 0; ++i)
            (*this).viewDiag(i) /= alpha;
    } else {
        for (IndexType i = numOffDiags(); i >= 0; --i)
            (*this).viewDiag(i) /= alpha;
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

    return _engine(row, col);
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

    return _engine(row, col);
}










// -- views --------------------------------------------------------------------

// hermitian views
template <typename FS>
const typename SbMatrix<FS>::ConstHermitianView
SbMatrix<FS>::hermitian() const
{
    return ConstHermitianView(_engine, _upLo);
}

template <typename FS>
typename SbMatrix<FS>::HermitianView
SbMatrix<FS>::hermitian()
{
    return HermitianView(_engine, _upLo);
}

// symmetric views
template <typename FS>
const typename SbMatrix<FS>::ConstView
SbMatrix<FS>::symmetric() const
{
    return ConstView(_engine, _upLo);
}

template <typename FS>
typename SbMatrix<FS>::View
SbMatrix<FS>::symmetric()
{
    return View(_engine, _upLo);
}

// triangular views
template <typename FS>
const typename SbMatrix<FS>::ConstTriangularView
SbMatrix<FS>::triangular() const
{
    return ConstTriangularView(_engine, _upLo);
}

template <typename FS>
typename SbMatrix<FS>::TriangularView
SbMatrix<FS>::triangular()
{
    return TriangularView(_engine, _upLo);
}

template <typename FS>
const typename SbMatrix<FS>::ConstVectorView
SbMatrix<FS>::viewDiag(IndexType diag) const
{
    if ((_upLo==Upper)^(diag>0))
        ConstVectorView(_engine.viewDiag(-diag, _engine.firstIndex()));
    
    return ConstVectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

template <typename FS>
typename SbMatrix<FS>::VectorView
SbMatrix<FS>::viewDiag(IndexType diag)
{
    if ((_upLo==Upper)^(diag>0))
         return VectorView(_engine.viewDiag(-diag, _engine.firstIndex()));
    
    return VectorView(_engine.viewDiag(diag, _engine.firstIndex()));
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::dim() const
{
    return _engine.dim();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numCols() const
{
    return _engine.numCols();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::numOffDiags() const
{
    return (_upLo==Upper) 
            ? _engine.numSuperDiags() 
            : _engine.numSubDiags();
}

template <typename FS>
typename SbMatrix<FS>::IndexType
SbMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
SbMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
const typename SbMatrix<FS>::ElementType *
SbMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename SbMatrix<FS>::ElementType *
SbMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
template <typename RHS>
bool
SbMatrix<FS>::resize(const SbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
SbMatrix<FS>::resize(IndexType dim, 
                     IndexType numOffDiags, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(dim, dim, 
                          (_upLo==Lower) ? numOffDiags : 0,  
                          (_upLo==Upper) ? numOffDiags : 0,  
                          firstIndex, value);
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename SbMatrix<FS>::Engine &
SbMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename SbMatrix<FS>::Engine &
SbMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
SbMatrix<FS>::upLo() const
{
    return _upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SBMATRIX_TCC
