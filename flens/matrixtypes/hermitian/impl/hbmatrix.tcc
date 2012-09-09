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

#include <flens/auxiliary/auxiliary.h>
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
    : _engine(dim, dim, (upLo == Lower) ? numOffDiags : 0, 
              (upLo == Upper) ? numOffDiags : 0, firstIndex),
      _upLo(upLo)
{
    ASSERT(dim>=0);
}

template <typename FS>
HbMatrix<FS>::HbMatrix(const Engine &engine, StorageUpLo upLo)
    : _engine(engine), _upLo(upLo)
{
}

template <typename FS>
HbMatrix<FS>::HbMatrix(const HbMatrix &rhs)
    : HermitianMatrix<HbMatrix<FS> >(),
      _engine(rhs.engine()),
      _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HbMatrix<FS>::HbMatrix(const HbMatrix<RHS> &rhs)
    : _engine(rhs.engine()),
      _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HbMatrix<FS>::HbMatrix(HbMatrix<RHS> &rhs)
    : _engine(rhs.engine()),
      _upLo(rhs.upLo())
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
HbMatrix<FS> &
HbMatrix<FS>::operator+=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HbMatrix<FS> &
HbMatrix<FS>::operator-=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
HbMatrix<FS> &
HbMatrix<FS>::operator*=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
     
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
HbMatrix<FS> &
HbMatrix<FS>::operator/=(const ElementType &alpha)
{
    ASSERT(cxxblas::imag(alpha)==0);
    
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
    return _engine(row, col);
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
    return _engine(row, col);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
typename HbMatrix<FS>::ConstView
HbMatrix<FS>::hermitian() const
{
    return ConstView(_engine, _upLo);
}

template <typename FS>
typename HbMatrix<FS>::View
HbMatrix<FS>::hermitian()
{
    return View(_engine, _upLo);
}

// symmetric view
template <typename FS>
typename HbMatrix<FS>::ConstSymmetricView
HbMatrix<FS>::symmetric() const
{
    return ConstSymmetricView(_engine, _upLo);
}

template <typename FS>
typename HbMatrix<FS>::SymmetricView
HbMatrix<FS>::symmetric()
{
    return SymmetricView(_engine, _upLo);
}

// triangular view
template <typename FS>
typename HbMatrix<FS>::ConstTriangularView
HbMatrix<FS>::triangular() const
{
    return ConstTriangularView(_engine, _upLo);
}

template <typename FS>
typename HbMatrix<FS>::TriangularView
HbMatrix<FS>::triangular()
{
    return TriangularView(_engine, _upLo);
}

template <typename FS>
const typename HbMatrix<FS>::ConstVectorView
HbMatrix<FS>::viewDiag(IndexType diag) const
{
    return ConstVectorView(_engine.viewDiag(diag));
}

template <typename FS>
typename HbMatrix<FS>::VectorView
HbMatrix<FS>::viewDiag(IndexType diag)
{
    return VectorView(_engine.viewDiag(diag));
}

// -- methods ------------------------------------------------------------------

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::dim() const
{
    return _engine.dim();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numCols() const
{
    return _engine.numCols();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numRows() const
{
    return _engine.numRows();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::leadingDimension() const
{
    return _engine.leadingDimension();
}

template <typename FS>
StorageOrder
HbMatrix<FS>::order() const
{
    return _engine.order;
}

template <typename FS>
typename HbMatrix<FS>::IndexType
HbMatrix<FS>::numOffDiags() const
{
    return (_upLo==Upper) ? _engine.numSuperDiags() 
                          : _engine.numSubDiags();
}

template <typename FS>
const typename HbMatrix<FS>::ElementType *
HbMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename HbMatrix<FS>::ElementType *
HbMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
template <typename RHS>
bool
HbMatrix<FS>::resize(const HbMatrix<RHS> &rhs,
                     const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
HbMatrix<FS>::resize(IndexType dim, IndexType numOffDiags,
                      IndexType firstIndex, const ElementType &value)
{
    return _engine.resize(dim, dim, 
                          (_upLo==Lower) ? numOffDiags : 0,  
                          (_upLo==Upper) ? numOffDiags : 0,  
                          firstIndex, value);
}

template <typename FS>
bool
HbMatrix<FS>::fill(const ElementType &value)
{
    ASSERT(cxxblas::imag(value)==0);
    
    return _engine.fill(value);
}

template <typename FS>
bool
HbMatrix<FS>::fillRandom()
{
    bool val = _engine.fillRandom();
    VectorView d = (*this).viewDiag();
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
    return _engine;
}

template <typename FS>
typename HbMatrix<FS>::Engine &
HbMatrix<FS>::engine()
{
    return _engine;
}

template <typename FS>
StorageUpLo
HbMatrix<FS>::upLo() const
{
    return _upLo;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HBMATRIX_TCC
