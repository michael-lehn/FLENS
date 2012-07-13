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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC 1

#include <flens/blas/level1/copy.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
HpMatrix<FS>::HpMatrix()
{
}

template <typename FS>
HpMatrix<FS>::HpMatrix(IndexType dim)
    : _engine(dim)
{
    ASSERT(dim>=0);
}

template <typename FS>
HpMatrix<FS>::HpMatrix(IndexType dim, IndexType firstIndex)
    : _engine(dim, dim, firstIndex)
{
    ASSERT(dim>=0);
}

template <typename FS>
HpMatrix<FS>::HpMatrix(const Engine &engine)
    : _engine(engine)
{
}

template <typename FS>
HpMatrix<FS>::HpMatrix(const HpMatrix &rhs)
    : HermitianMatrix<HpMatrix<FS> >(),
      _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
HpMatrix<FS>::HpMatrix(const HpMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
HpMatrix<FS>::HpMatrix(HpMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
HpMatrix<FS>::HpMatrix(const Matrix<RHS> &rhs)
{
    blas::copy(rhs.impl(), *this);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
HpMatrix<FS> &
HpMatrix<FS>::operator=(const HpMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
HpMatrix<FS> &
HpMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
HpMatrix<FS> &
HpMatrix<FS>::operator=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HpMatrix<FS> &
HpMatrix<FS>::operator+=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HpMatrix<FS> &
HpMatrix<FS>::operator-=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HpMatrix<FS> &
HpMatrix<FS>::operator*=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HpMatrix<FS> &
HpMatrix<FS>::operator/=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
const typename HpMatrix<FS>::ElementType &
HpMatrix<FS>::operator()(IndexType row, IndexType col) const
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
typename HpMatrix<FS>::ElementType &
HpMatrix<FS>::operator()(IndexType row, IndexType col)
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

// general views
template <typename FS>
typename HpMatrix<FS>::ConstView
HpMatrix<FS>::hermitian() const
{
    return ConstView(_engine);
}

template <typename FS>
typename HpMatrix<FS>::View
HpMatrix<FS>::hermitian()
{
    return View(_engine);
}

// symmetric view
template <typename FS>
typename HpMatrix<FS>::ConstSymmetricView
HpMatrix<FS>::symmetric() const
{
    return ConstSymmetricView(_engine);
}

template <typename FS>
typename HpMatrix<FS>::SymmetricView
HpMatrix<FS>::symmetric()
{
    return SymmetricView(_engine);
}

// triangular view
template <typename FS>
typename HpMatrix<FS>::ConstTriangularView
HpMatrix<FS>::triangular() const
{
    return ConstTriangularView(_engine);
}

template <typename FS>
typename HpMatrix<FS>::TriangularView
HpMatrix<FS>::triangular()
{
    return TriangularView(_engine);
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename HpMatrix<FS>::IndexType
HpMatrix<FS>::dim() const
{
    return _engine.dim();
}

template <typename FS>
typename HpMatrix<FS>::IndexType
HpMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename HpMatrix<FS>::IndexType
HpMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
const typename HpMatrix<FS>::ElementType *
HpMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename HpMatrix<FS>::ElementType *
HpMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
template <typename RHS>
bool
HpMatrix<FS>::resize(const HpMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
HpMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                     const ElementType &value)
{
    return _engine.resize(dim, firstIndex, value);
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const typename HpMatrix<FS>::Engine &
HpMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename HpMatrix<FS>::Engine &
HpMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
HpMatrix<FS>::upLo() const
{
    return FS::upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_TCC
