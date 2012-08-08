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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIXINITIALIZER_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIXINITIALIZER_TCC 1

#include <flens/auxiliary/macros.h>
#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/general/impl/gematrixinitializer.h>

namespace flens {

template <typename FS>
GeMatrixInitializer<FS>::GeMatrixInitializer(GeMatrix<FS> &A,
                                             IndexType row,
                                             IndexType col)
    : _A(A), _row(row), _col(col)
{
}

template <typename FS>
GeMatrixInitializer<FS>
GeMatrixInitializer<FS>::operator,(const ElementType &value)
{
    ++_col;
    if (_col>_A.lastCol()) {
        _col = _A.firstCol();
        ++_row;
    }
    ASSERT(_col>=_A.firstCol());
    ASSERT(_col<=_A.lastCol());
    ASSERT(_row>=_A.firstRow());
    ASSERT(_row<=_A.lastRow());
    _A(_row, _col) = value;
    return GeMatrixInitializer(_A, _row, _col);
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIXINITIALIZER_TCC
