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

#ifndef FLENS_VECTORTYPES_IMPL_TINYVECTOR_TCC
#define FLENS_VECTORTYPES_IMPL_TINYVECTOR_TCC 1

#include <flens/blas/blas.h>
#include <flens/typedefs.h>

#include <flens/vectortypes/impl/tinyvector.h>

namespace flens {

// -- constructors -------------------------------------------------------------

template <typename TA>
TinyVector<TA>::TinyVector()
{
}

template <typename TA>
TinyVector<TA>::TinyVector(const Engine &engine)
    : _array(engine)
{
}

template <typename TA>
template <typename RHS>
TinyVector<TA>::TinyVector(const TinyVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename TA>
template <typename RHS>
TinyVector<TA>::TinyVector(TinyVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename TA>
template <typename RHS>
TinyVector<TA>::TinyVector(const Vector<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- operators -----------------------------------------------------------------

template <typename TA>
TinyVector<TA> &
TinyVector<TA>::operator=(const TinyVector &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename TA>
template <typename RHS>
TinyVector<TA> &
TinyVector<TA>::operator=(const Vector<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename TA>
template <typename RHS>
TinyVector<TA> &
TinyVector<TA>::operator+=(const Vector<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename TA>
template <typename RHS>
TinyVector<TA> &
TinyVector<TA>::operator-=(const Vector<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename TA>
TinyVector<TA> &
TinyVector<TA>::operator+=(const ElementType &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) += alpha;
    }
    return *this;
}

template <typename TA>
TinyVector<TA> &
TinyVector<TA>::operator-=(const ElementType &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) -= alpha;
    }
    return *this;
}

template <typename TA>
TinyVector<TA> &
TinyVector<TA>::operator*=(const ElementType &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename TA>
TinyVector<TA> &
TinyVector<TA>::operator/=(const ElementType &alpha)
{
    blas::rscal(alpha, *this);
    return *this;
}

template <typename TA>
const typename TinyVector<TA>::ElementType &
TinyVector<TA>::operator()(IndexType index) const
{
    return _array(index);
}

template <typename TA>
typename TinyVector<TA>::ElementType &
TinyVector<TA>::operator()(IndexType index)
{
    return _array(index);
}

//-- methods -------------------------------------------------------------------

template <typename TA>
typename TinyVector<TA>::IndexType
TinyVector<TA>::length() const
{
    return Engine::length;
}

template <typename TA>
typename TinyVector<TA>::IndexType
TinyVector<TA>::firstIndex() const
{
    return Engine::firstIndex;
}

template <typename TA>
typename TinyVector<TA>::IndexType
TinyVector<TA>::lastIndex() const
{
    return Engine::lastIndex;
}

template <typename TA>
const typename TinyVector<TA>::ElementType *
TinyVector<TA>::data() const
{
    return _array.data();
}

template <typename TA>
typename TinyVector<TA>::ElementType *
TinyVector<TA>::data()
{
    return _array.data();
}

template <typename TA>
typename TinyVector<TA>::IndexType
TinyVector<TA>::stride() const
{
    return Engine::stride;
}

template <typename TA>
void
TinyVector<TA>::fill(const ElementType &value)
{
    engine().fill(value);
}

//-- implementation ------------------------------------------------------------
template <typename TA>
const typename TinyVector<TA>::Engine &
TinyVector<TA>::engine() const
{
    return _array;
}

template <typename TA>
typename TinyVector<TA>::Engine &
TinyVector<TA>::engine()
{
    return _array;
}

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_TINYVECTOR_TCC
