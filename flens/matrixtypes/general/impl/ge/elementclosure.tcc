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

#include <flens/matrixtypes/general/impl/ge/elementclosure.h>
#include <flens/typedefs.h>

namespace flens { namespace gematrix {

template <typename M>
ElementClosure<M>::ElementClosure(Matrix &matrix,
                                  IndexVariable &row, IndexVariable &col)
    : matrix_(matrix), row_(row), col_(col)
{
}

template <typename M>
//void
int
ElementClosure<M>::operator=(const ElementType &rhs)
{
    typename IndexVariable::ElementType &i = row_.value();
    typename IndexVariable::ElementType &j = col_.value();

    if (M::Engine::order==RowMajor) {
        for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
            for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
                value() = rhs;
            }
        }
    } else {
        for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
            for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
                value() = rhs;
            }
        }
    }
    return 0;
}

template <typename M>
template <typename S>
void
ElementClosure<M>::operator=(const Scalar<S> &rhs)
{
    typename IndexVariable::ElementType &i = row_.value();
    typename IndexVariable::ElementType &j = col_.value();

    if (M::Engine::order==RowMajor) {
        for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
            for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
                value() = rhs.impl().value();
            }
        }
    } else {
        for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
            for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
                value() = rhs.impl().value();
            }
        }
    }
}

template <typename M>
void
ElementClosure<M>::operator=(const ElementClosure &rhs)
{
    typename IndexVariable::ElementType &i = row_.value();
    typename IndexVariable::ElementType &j = col_.value();

    if (M::Engine::order==RowMajor) {
        for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
            for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
                value() = rhs.impl().value();
            }
        }
    } else {
        for (j=matrix_.firstCol(); j<=matrix_.lastCol(); ++j) {
            for (i=matrix_.firstRow(); i<=matrix_.lastRow(); ++i) {
                value() = rhs.impl().value();
            }
        }
    }
}

template <typename M>
const typename ElementClosure<M>::ElementType &
ElementClosure<M>::value() const
{
    return matrix_(row_.value(), col_.value());
}

template <typename M>
typename ElementClosure<M>::ElementType &
ElementClosure<M>::value()
{
    return matrix_(row_.value(), col_.value());
}

} } // namespace gematrix, flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_TCC
