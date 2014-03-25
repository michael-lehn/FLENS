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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/level1/copy.h>
#include <flens/matrixtypes/general/impl/diagmatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename FS>
DiagMatrix<FS>::DiagMatrix()
{
}

template <typename FS>
DiagMatrix<FS>::DiagMatrix(IndexType dim, IndexType firstIndex)
      : _engine(dim, firstIndex)
{
}

template <typename FS>
DiagMatrix<FS>::DiagMatrix(const Engine &engine)
    : _engine(engine)
{
}

template <typename FS>
DiagMatrix<FS>::DiagMatrix(const DiagMatrix &rhs)
    : GeneralMatrix<DiagMatrix<FS> >(), _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS>::DiagMatrix(const DiagMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS>::DiagMatrix(DiagMatrix<RHS> &rhs)
    : _engine(rhs.engine())
{
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS>::DiagMatrix(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS>::DiagMatrix(const DenseVector<RHS> &rhs)
    : _engine(rhs.engine())
{

}

// -- operators ----------------------------------------------------------------

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator=(const DiagMatrix &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator=(const Matrix<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator=(const DenseVector<RHS> &rhs)
{
    (*this).diag() = rhs;
    return *this;
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator+=(const Matrix<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename FS>
template <typename RHS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator-=(const Matrix<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator=(const ElementType &alpha)
{
    (*this).diag() = alpha;
    
    return *this;
}

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator+=(const ElementType &alpha)
{
    (*this).diag() += alpha;
    
    return *this;
}

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator-=(const ElementType &alpha)
{
    
    (*this).diag() -= alpha;
    
    return *this;
}

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator*=(const ElementType &alpha)
{
   (*this).diag() *= alpha;
    
    return *this;
}

template <typename FS>
DiagMatrix<FS> &
DiagMatrix<FS>::operator/=(const ElementType &alpha)
{
    (*this).diag() /= alpha;
    
    return *this;
}

template <typename FS>
const typename DiagMatrix<FS>::ElementType &
DiagMatrix<FS>::operator()(IndexType row, IndexType DEBUG_VAR(col)) const
{
    ASSERT( row==col );

    return _engine(row);
}

template <typename FS>
typename DiagMatrix<FS>::ElementType &
DiagMatrix<FS>::operator()(IndexType row, IndexType DEBUG_VAR(col))
{
    ASSERT( row==col );
    
    return _engine(row);
}

// -- views --------------------------------------------------------------------
template <typename FS>
const typename DiagMatrix<FS>::ConstVectorView
DiagMatrix<FS>::diag() const
{

    return ConstVectorView(_engine.view(_engine.firstIndex(), 
                                        _engine.lastIndex(), 
                                        _engine.stride()));
}

template <typename FS>
typename DiagMatrix<FS>::VectorView
DiagMatrix<FS>::diag()
{
    return VectorView(_engine.view(_engine.firstIndex(),
                                   _engine.lastIndex(),
                                   _engine.stride()));
}

// -- methods ------------------------------------------------------------------
template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::dim() const
{
    return _engine.length();
}

template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::numCols() const
{
    return _engine.length();
}

template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::numRows() const
{
    return _engine.length();
}

template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::firstIndex() const
{
    return _engine.firstIndex();
}

template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::lastIndex() const
{
    return _engine.lastIndex();
}

template <typename FS>
typename DiagMatrix<FS>::IndexType
DiagMatrix<FS>::leadingDimension() const
{
    return _engine.stride();
}

template <typename FS>
const typename DiagMatrix<FS>::ElementType *
DiagMatrix<FS>::data() const
{
    return _engine.data();
}

template <typename FS>
typename DiagMatrix<FS>::ElementType *
DiagMatrix<FS>::data()
{
    return _engine.data();
}

template <typename FS>
template <typename RHS>
bool
DiagMatrix<FS>::resize(const DiagMatrix<RHS> &rhs,
                       const ElementType &value)
{
    return _engine.resize(rhs.engine(), value);
}

template <typename FS>
bool
DiagMatrix<FS>::resize(IndexType dim, IndexType firstIndex,
                       const ElementType &value)
{
    return _engine.resize(dim, firstIndex, value);
}

template <typename FS>
bool
DiagMatrix<FS>::fill(const ElementType &value)
{
    return _engine.fill(value);
}


// -- implementation -----------------------------------------------------------

template <typename FS>
const typename DiagMatrix<FS>::Engine &
DiagMatrix<FS>::engine() const
{
    return _engine;
}

template <typename FS>
typename DiagMatrix<FS>::Engine &
DiagMatrix<FS>::engine()
{
    return _engine;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_TCC