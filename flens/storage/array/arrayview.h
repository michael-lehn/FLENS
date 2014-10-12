/*
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

#ifndef FLENS_STORAGE_ARRAY_ARRAYVIEW_H
#define FLENS_STORAGE_ARRAY_ARRAYVIEW_H 1

#include <cxxstd/memory.h>
#include <flens/storage/indexoptions.h>

namespace flens {

template <typename T, typename I, typename A>
    class ConstArrayView;

template <typename T, typename I, typename A>
    class Array;

template <typename T,
          typename I = IndexOptions<>,
          typename A = std::allocator<T> >
class ArrayView
{
    public:
        typedef T                       ElementType;
        typedef typename I::IndexType   IndexType;
        typedef A                       Allocator;

        typedef ConstArrayView<T, I, A> ConstView;
        typedef ArrayView               View;
        typedef Array<T, I, A>          NoView;

        static const IndexType          defaultIndexBase = I::defaultIndexBase;

        ArrayView(IndexType length,
                  ElementType *data,
                  IndexType stride = IndexType(1),
                  IndexType firstIndex = defaultIndexBase,
                  const Allocator &allocator = Allocator());

        ArrayView(const ArrayView &rhs);

        template <typename RHS>
            ArrayView(RHS &rhs);

        ~ArrayView();

        //-- operators ---------------------------------------------------------

        const ElementType &
        operator()(IndexType index) const;

        ElementType &
        operator()(IndexType index);

        //-- methods -----------------------------------------------------------

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        IndexType
        length() const;

        IndexType
        stride() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        const Allocator &
        allocator() const;

        bool
        resize(IndexType length,
               IndexType firstIndex = defaultIndexBase,
               const ElementType &value = ElementType());

        template <typename ARRAY>
            bool
            resize(const ARRAY &rhs, const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        void
        changeIndexBase(IndexType firstIndex);

        const ConstView
        view(IndexType from, IndexType to,
             IndexType stride = IndexType(1),
             IndexType firstViewIndex =  defaultIndexBase) const;

        ArrayView
        view(IndexType from, IndexType to,
             IndexType stride = IndexType(1),
             IndexType firstViewIndex = defaultIndexBase);

    private:
        ElementType  *data_;
        Allocator    allocator_;
        IndexType    length_, stride_, firstIndex_;
};

//-- ArrayView specific functions ----------------------------------------------

//
//  fillRandom
//

template <typename T, typename I, typename A>
    bool
    fillRandom(ArrayView<T, I, A> &x);

} // namespace flens

#endif // FLENS_STORAGE_ARRAY_ARRAYVIEW_H
