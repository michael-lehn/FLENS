/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_SCALARCLOSURES_H
#define FLENS_SCALARCLOSURES_H 1

#include <flens/matvec.h>
#include <flens/operationtypes.h>
#include <flens/traits.h>

namespace flens {

// == ScalarClosure ============================================================

template <typename E>
class ScalarClosure
{
    public:
        typedef E                                Engine;
        typedef typename Engine::ElementType     ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        ScalarClosure();

        ScalarClosure(const Engine &engine);

        virtual
        ~ScalarClosure();

        operator ElementType() const
        {
            return _engine.value();
        }

        void
        operator=(const ScalarClosure<Engine> &rhs);

        template <typename RHS>
            void
            operator=(const ScalarClosure<RHS> &rhs);

        const Engine &
        engine() const;

        Engine &
        engine();

    private:
        Engine _engine;
};

//== SimpleIndexEngine =========================================================

class SimpleIndexEngine
{
    public:
        typedef int                     ElementType;
        typedef Ref<ElementType>::Type  ElementTypeRef;

        void
        setRange(int firstIndex, int lastIndex);

        int
        firstIndex() const;

        int
        lastIndex() const;

        int
        value() const;

        int &
        value();

    private:
        int _firstIndex;
        int _lastIndex;
        int _value;
};

typedef ScalarClosure<SimpleIndexEngine> Index;

//== Scalar: Wrapper for scalar values =========================================

template <typename T>
class Scalar
{
    public:
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        Scalar(T value);

        operator T() const;

        const T&
        value() const;

    private:
        T _value;
};

template <typename T>
struct Ref<Scalar<T> >
{
    typedef typename Ref<T>::Type Type;
};

// == ScalarOperation ==========================================================

template <typename Op, typename L, typename R>
class ScalarOperation
{
    public:
        typedef typename L::ElementType          TL;
        typedef typename R::ElementType          TR;
        typedef typename Promotion<TL,TR>::Type  T;
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        ScalarOperation(const L &l, const R &r);

        T
        value() const;

    private:
        const L &_left;
        const R &_right;
};

// == ScalarUnaryOperation =====================================================

template <typename Op, typename Arg>
class ScalarUnaryOperation
{
    public:
        typedef typename Arg::ElementType        T;
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        ScalarUnaryOperation(const Arg &arg);

        T
        value() const;

    private:
        const Arg &_arg;
};

// == VectorElement ============================================================

template <typename VecImpl>
class VectorElement
{
    public:
        typedef typename VecImpl::ElementType    T;
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        VectorElement(VecImpl &vecImpl, Index &index);

        const T &
        value () const;

        VecImpl &_vecImpl;
        Index &_index;
};

//------------------------------------------------------------------------------

template <typename X, typename Y>
    void
    copy(const ScalarClosure<X> &x, VectorElement<Y> &y);

// == MatrixElement ============================================================

template <typename MatImpl>
class MatrixElement
{
    public:
        typedef typename MatImpl::ElementType    T;
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;


        MatrixElement(MatImpl &matImpl, Index &i, Index &j);

        const T &
        value () const;

        MatImpl &_matImpl;
        Index   &_i;
        Index   &_j;
};

//------------------------------------------------------------------------------

template <typename X, typename Y>
void
copy(const ScalarClosure<X> &x, MatrixElement<Y> &y);

// == Summation ================================================================

template <typename I>
class Summation
{
    public:
        typedef typename I::ElementType          T;
        typedef T                                ElementType;
        typedef typename Ref<ElementType>::Type  ElementTypeRef;

        Summation(const ScalarClosure<I> &summand, Index &index);

        T
        value() const;

    private:
        const I &_summand;
        Index &_index;
};

//------------------------------------------------------------------------------

template <typename I>
    const ScalarClosure<Summation<I> >
    sum(const ScalarClosure<I> &summand, Index &index);

//==============================================================================

} // namespace flens

#include <flens/scalarclosures.tcc>

#endif // FLENS_SCALARCLOSURES_H
