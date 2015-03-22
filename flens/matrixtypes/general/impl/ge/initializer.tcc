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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GE_INITIALIZER_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GE_INITIALIZER_TCC 1

#include <flens/auxiliary/macros.h>
#include <flens/matrixtypes/general/impl/gematrix.h>

namespace flens { namespace gematrix {

template <typename M>
Initializer<M>::Initializer(Matrix &A, IndexType row, IndexType col)
    : A_(A), row_(row), col_(col)
{
}

template <typename M>
Initializer<M>
Initializer<M>::operator,(const T &value)
{
    ++col_;
    if (col_>A_.lastCol()) {
        col_ = A_.firstCol();
        ++row_;
    }
    ASSERT(col_>=A_.firstCol());
    ASSERT(col_<=A_.lastCol());
    ASSERT(row_>=A_.firstRow());
    ASSERT(row_<=A_.lastRow());
    A_(row_, col_) = value;
    return Initializer(A_, row_, col_);
}

} } // namespace gematrix, flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GE_INITIALIZER_TCC
