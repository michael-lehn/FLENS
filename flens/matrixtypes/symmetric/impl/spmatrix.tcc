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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_TCC
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_TCC 1

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/symmetric/impl/symatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
SpMatrix<FS>::SpMatrix()
{
  
}
  
template <typename FS>
SpMatrix<FS>::SpMatrix(IndexType dim)
      : _engine(dim)
{
  
}
  
template <typename FS>
SpMatrix<FS>::SpMatrix(const Engine &engine)
    : _engine(engine)
{
}

template <typename FS>
SpMatrix<FS>::SpMatrix(const SpMatrix &rhs)
    : SymmetricMatrix<SpMatrix<FS> >(),
      _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
SpMatrix<FS>::SpMatrix(const SpMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
SpMatrix<FS>::SpMatrix(SpMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
SpMatrix<FS>::SpMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator=(const SpMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
SpMatrix<FS> &
SpMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator=(const ElementType &alpha)
{
    if (upLo()==Upper) {
        for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = i; j <= lastIndex(); ++j)
                (*this)(i,j) = alpha;
 
    } else {
      for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = firstIndex(); j <= i; ++j)
                (*this)(i,j) = alpha;
    }
    return *this;
}

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator+=(const ElementType &alpha)
{

    if (upLo()==Upper) {
        for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = i; j <= lastIndex(); ++j)
                (*this)(i,j) += alpha;
 
    } else {
      for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = firstIndex(); j <= i; ++j)
                (*this)(i,j) += alpha;
    }
    return *this;
}

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator-=(const ElementType &alpha)
{

    if (upLo()==Upper) {
        for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = i; j <= lastIndex(); ++j)
                (*this)(i,j) -= alpha;
 
    } else {
      for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = firstIndex(); j <= i; ++j)
                (*this)(i,j) -= alpha;
    }
    return *this;
}

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator*=(const ElementType &alpha)
{
   
    if (upLo()==Upper) {
        for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = i; j <= lastIndex(); ++j)
                (*this)(i,j) *= alpha;
 
    } else {
      for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = firstIndex(); j <= i; ++j)
                (*this)(i,j) *= alpha;
    }
    return *this;
}

template <typename FS>
SpMatrix<FS> &
SpMatrix<FS>::operator/=(const ElementType &alpha)
{
   
    if (upLo()==Upper) {
        for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = i; j <= lastIndex(); ++j)
                (*this)(i,j) /= alpha;
 
    } else {
      for (IndexType i = firstIndex(); i<= lastIndex(); ++i)
            for(IndexType j = firstIndex(); j <= i; ++j)
                (*this)(i,j) /= alpha;
    }
    return *this;
}


template <typename FS>
const typename SpMatrix<FS>::ElementType &
SpMatrix<FS>::operator()(IndexType row, IndexType col) const
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
typename SpMatrix<FS>::ElementType &
SpMatrix<FS>::operator()(IndexType row, IndexType col)
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
const typename SpMatrix<FS>::ConstHermitianView
SpMatrix<FS>::hermitian() const
{
    return ConstHermitianView(_engine);
}

template <typename FS>
typename SpMatrix<FS>::HermitianView
SpMatrix<FS>::hermitian()
{
    return HermitianView(_engine);
}

// symmetric views
template <typename FS>
const typename SpMatrix<FS>::ConstView
SpMatrix<FS>::symmetric() const
{
    return ConstView(_engine);
}

template <typename FS>
typename SpMatrix<FS>::View
SpMatrix<FS>::symmetric()
{
    return View(_engine);
}

// triangular views
template <typename FS>
const typename SpMatrix<FS>::ConstTriangularView
SpMatrix<FS>::triangular() const
{
    return ConstTriangularView(_engine);
}

template <typename FS>
typename SpMatrix<FS>::TriangularView
SpMatrix<FS>::triangular()
{
    return TriangularView(_engine);
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename SpMatrix<FS>::IndexType
SpMatrix<FS>::dim() const
{
    return _engine.dim();
}

template <typename FS>
typename SpMatrix<FS>::IndexType
SpMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename SpMatrix<FS>::IndexType
SpMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
const typename SpMatrix<FS>::ElementType *
SpMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename SpMatrix<FS>::ElementType *
SpMatrix<FS>::data()
{
    return _engine.data();
}


template <typename FS>
StorageOrder
SpMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
template <typename RHS>
bool
SpMatrix<FS>::resize(const SpMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
SpMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(dim, firstIndex, value);
}

template <typename FS>
bool
SpMatrix<FS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}

template <typename FS>
bool
SpMatrix<FS>::fillRandom()
{
    return _engine.fillRandom();
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename SpMatrix<FS>::Engine &
SpMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename SpMatrix<FS>::Engine &
SpMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
SpMatrix<FS>::upLo() const
{
    return FS::upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_TCC
