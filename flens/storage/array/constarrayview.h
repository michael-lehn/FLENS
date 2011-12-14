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

#ifndef FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_H
#define FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_H 1

namespace flens {

template <typename T, typename I, typename A>
    class ArrayView;

template <typename T, typename I, typename A>
    class Array;

template <typename T,
          typename I = IndexOptions<>,
          typename A = std::allocator<T> >
class ConstArrayView
{
    public:
        typedef T                      ElementType;
        typedef typename I::IndexType  IndexType;
        typedef A                      Allocator;

        typedef ConstArrayView          ConstView;
        typedef ArrayView<T, I, A>      View;
        typedef Array<T, I, A>          NoView;

        static const IndexType          defaultIndexBase = I::defaultIndexBase;

        ConstArrayView(IndexType length,
                       const ElementType *data,
                       IndexType stride = IndexType(1),
                       IndexType firstIndex =  defaultIndexBase,
                       const Allocator &allocator = Allocator());

        ConstArrayView(const ConstArrayView &rhs);

        template <typename RHS>
            ConstArrayView(const RHS &rhs);

        ~ConstArrayView();

        //-- operators ---------------------------------------------------------
        const ElementType &
        operator()(IndexType index) const;

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

        const Allocator &
        allocator() const;

        const ConstArrayView
        view(IndexType from, IndexType to,
             IndexType stride = IndexType(1),
             IndexType firstViewIndex =  defaultIndexBase) const;

    private:
        ConstArrayView &
        operator=(const ConstArrayView &rhs);

        const ElementType *_data;
        Allocator         _allocator;
        IndexType         _length, _stride, _firstIndex;
};

} // namespace flens

#endif // FLENS_STORAGE_ARRAY_CONSTARRAYVIEW_H
