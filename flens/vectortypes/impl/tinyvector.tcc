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

template <typename A>
TinyVector<A>::TinyVector()
{
}

template <typename A>
TinyVector<A>::TinyVector(const Engine &engine)
    : _array(engine)
{
}

template <typename A>
template <typename RHS>
TinyVector<A>::TinyVector(const TinyVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename A>
template <typename RHS>
TinyVector<A>::TinyVector(TinyVector<RHS> &rhs)
    : _array(rhs.engine())
{
}

template <typename A>
template <typename RHS>
TinyVector<A>::TinyVector(const Vector<RHS> &rhs)
{
    assign(rhs, *this);
}

//-- operators -----------------------------------------------------------------

template <typename A>
TinyVector<A> &
TinyVector<A>::operator=(const TinyVector &rhs)
{
    if (this!=&rhs) {
        assign(rhs, *this);
    }
    return *this;
}

template <typename A>
template <typename RHS>
TinyVector<A> &
TinyVector<A>::operator=(const Vector<RHS> &rhs)
{
    assign(rhs, *this);
    return *this;
}

template <typename A>
template <typename RHS>
TinyVector<A> &
TinyVector<A>::operator+=(const Vector<RHS> &rhs)
{
    plusAssign(rhs, *this);
    return *this;
}

template <typename A>
template <typename RHS>
TinyVector<A> &
TinyVector<A>::operator-=(const Vector<RHS> &rhs)
{
    minusAssign(rhs, *this);
    return *this;
}

template <typename A>
TinyVector<A> &
TinyVector<A>::operator+=(const ElementType &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) += alpha;
    }
    return *this;
}

template <typename A>
TinyVector<A> &
TinyVector<A>::operator-=(const ElementType &alpha)
{
    for (int i=firstIndex(); i<=lastIndex(); ++i) {
        (*this)(i) -= alpha;
    }
    return *this;
}

template <typename A>
TinyVector<A> &
TinyVector<A>::operator*=(const ElementType &alpha)
{
    blas::scal(alpha, *this);
    return *this;
}

template <typename A>
TinyVector<A> &
TinyVector<A>::operator/=(const ElementType &alpha)
{
    blas::rscal(alpha, *this);
    return *this;
}

template <typename A>
const typename TinyVector<A>::ElementType &
TinyVector<A>::operator()(IndexType index) const
{
    return _array(index);
}

template <typename A>
typename TinyVector<A>::ElementType &
TinyVector<A>::operator()(IndexType index)
{
    return _array(index);
}

//-- methods -------------------------------------------------------------------

template <typename A>
typename TinyVector<A>::IndexType
TinyVector<A>::firstIndex() const
{
    return Engine::firstIndex;
}

template <typename A>
typename TinyVector<A>::IndexType
TinyVector<A>::lastIndex() const
{
    return Engine::lastIndex;
}

template <typename A>
typename TinyVector<A>::IndexType
TinyVector<A>::length() const
{
    return Engine::length;
}

template <typename A>
const typename TinyVector<A>::ElementType *
TinyVector<A>::data() const
{
    return _array.data();
}

template <typename A>
typename TinyVector<A>::ElementType *
TinyVector<A>::data()
{
    return _array.data();
}

template <typename A>
typename TinyVector<A>::IndexType
TinyVector<A>::stride() const
{
    return Engine::stride;
}

//-- implementation ------------------------------------------------------------
template <typename A>
const typename TinyVector<A>::Engine &
TinyVector<A>::engine() const
{
    return _array;
}

template <typename A>
typename TinyVector<A>::Engine &
TinyVector<A>::engine()
{
    return _array;
}

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_TINYVECTOR_TCC
