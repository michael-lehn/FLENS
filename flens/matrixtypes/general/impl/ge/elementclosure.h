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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_H
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_H 1

#include <flens/auxiliary/range.h>
#include <flens/scalartypes/scalar.h>
#include <flens/vectortypes/vector.h>

namespace flens { namespace gematrix {

template <typename M>
class ElementClosure
    : public Scalar<ElementClosure<M> >
{
    public:
        typedef M                               Matrix;
        typedef typename Matrix::ElementType    ElementType;
        typedef typename Matrix::IndexVariable  IndexVariable;

        ElementClosure(Matrix &matrix, IndexVariable &row, IndexVariable &col);

        //void
        int
        operator=(const ElementType &rhs);

        template <typename S>
            void
            operator=(const Scalar<S> &rhs);

        void
        operator=(const ElementClosure &rhs);

        const ElementType &
        value() const;

        ElementType &
        value();

    private:
        Matrix         &_matrix;
        IndexVariable  &_row;
        IndexVariable  &_col;
};

} } // namespace gematrix, flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GE_ELEMENTCLOSURE_H
