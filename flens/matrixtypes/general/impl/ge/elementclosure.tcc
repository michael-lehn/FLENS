/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_TCC
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_TCC 1

#include <cxxblas/cxxblas.h>
#include <flens/storage/storage.h>

namespace flens { namespace gematrix {

template <typename M>
ElementClosure<M>::ElementClosure(Matrix &matrix,
                                  IndexVariable &row, IndexVariable &col)
    : _matrix(matrix), _row(row), _col(col)
{
}

template <typename M>
void
ElementClosure<M>::operator=(const ElementType &rhs)
{
    typename IndexVariable::ElementType &i = _row.value();
    typename IndexVariable::ElementType &j = _col.value();

    if (StorageInfo<typename M::Engine>::Order==cxxblas::RowMajor) {
        for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
            for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
                value() = rhs;
            }
        }
    } else {
        for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
            for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
                value() = rhs;
            }
        }
    }
}

template <typename M>
template <typename S>
void
ElementClosure<M>::operator=(const Scalar<S> &rhs)
{
    typename IndexVariable::ElementType &i = _row.value();
    typename IndexVariable::ElementType &j = _col.value();

    if (StorageInfo<typename M::Engine>::Order==cxxblas::RowMajor) {
        for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
            for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
                value() = rhs.impl().value();
            }
        }
    } else {
        for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
            for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
                value() = rhs.impl().value();
            }
        }
    }
}

template <typename M>
void
ElementClosure<M>::operator=(const ElementClosure &rhs)
{
    typename IndexVariable::ElementType &i = _row.value();
    typename IndexVariable::ElementType &j = _col.value();

    if (StorageInfo<typename M::Engine>::Order==cxxblas::RowMajor) {
        for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
            for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
                value() = rhs.impl().value();
            }
        }
    } else {
        for (j=_matrix.firstCol(); j<=_matrix.lastCol(); ++j) {
            for (i=_matrix.firstRow(); i<=_matrix.lastRow(); ++i) {
                value() = rhs.impl().value();
            }
        }
    }
}

template <typename M>
const typename ElementClosure<M>::ElementType &
ElementClosure<M>::value() const
{
    return _matrix(_row.value(), _col.value());
}

template <typename M>
typename ElementClosure<M>::ElementType &
ElementClosure<M>::value()
{
    return _matrix(_row.value(), _col.value());
}

} } // namespace gematrix, flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_TCC
