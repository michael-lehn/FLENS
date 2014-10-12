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

#ifndef FLENS_VECTORTYPES_IMPL_DV_INITIALIZER_TCC
#define FLENS_VECTORTYPES_IMPL_DV_INITIALIZER_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/auxiliary/macros.h>
#include <flens/vectortypes/impl/densevector.h>

namespace flens { namespace densevector {

template <typename V>
Initializer<V>::Initializer(Vector &x, IndexType index)
    : x_(x), index_(index)
{
}

template <typename A>
Initializer<A>
Initializer<A>::operator,(const ElementType &value)
{
    index_ += x_.inc();
    ASSERT(index_>=std::min(x_.firstIndex(), x_.lastIndex()));
    ASSERT(index_<=std::max(x_.firstIndex(), x_.lastIndex()));
    x_(index_) = value;
    return Initializer(x_, index_);
}

} } // namespace densevector, flens

#endif // FLENS_VECTORTYPES_IMPL_DV_INITIALIZER_TCC
