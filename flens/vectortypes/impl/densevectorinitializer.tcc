/*
 *   Copyright (c) 2003, Alexander Stippler
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_VECTORTYPES_IMPL_DENSEVECTORINITIALIZER_TCC
#define FLENS_VECTORTYPES_IMPL_DENSEVECTORINITIALIZER_TCC 1

#include <flens/auxiliary/macros.h>
#include <flens/vectortypes/impl/densevector.h>
#include <flens/vectortypes/impl/densevectorinitializer.h>

namespace flens {

template <typename A>
DenseVectorInitializer<A>::DenseVectorInitializer(DenseVector<A> &x,
                                                  IndexType index)
    : _x(x), _index(index)
{
}

template <typename A>
DenseVectorInitializer<A>
DenseVectorInitializer<A>::operator,(const ElementType &value)
{
    ++_index;
    ASSERT(_index>=_x.firstIndex());
    ASSERT(_index<=_x.lastIndex());
    _x(_index) = value;
    return DenseVectorInitializer(_x, _index);
}

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_DENSEVECTORINITIALIZER_H
