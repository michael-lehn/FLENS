/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef FLENS_VECTORTYPES_IMPL_DENSEVECTOR_H
#define FLENS_VECTORTYPES_IMPL_DENSEVECTOR_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/scalartypes/scalar.h>
#include <flens/vectortypes/auxiliary/imagvector.h>
#include <flens/vectortypes/auxiliary/realvector.h>
#include <flens/storage/array/imagarray.h>
#include <flens/storage/array/realarray.h>
#include <flens/vectortypes/impl/dv/constelementclosure.h>
#include <flens/vectortypes/impl/dv/elementclosure.h>
#include <flens/vectortypes/impl/dv/initializer.h>
#include <flens/vectortypes/vector.h>

namespace flens {

// forward declaration
template <typename IndexType>
    class IndexVariable;

template <typename A>
class DenseVector
    : public Vector<DenseVector<A> >
{
    public:
        typedef A                                           Engine;
        typedef typename Engine::ElementType                ElementType;
        typedef typename Engine::IndexType                  IndexType;

        // view types from Engine
        typedef typename Engine::ConstView                  EngineConstView;
        typedef typename Engine::View                       EngineView;
        typedef typename Engine::NoView                     EngineNoView;

        // view types
        typedef DenseVector<EngineConstView>                ConstView;
        typedef DenseVector<EngineView>                     View;
        typedef DenseVector<EngineNoView>                   NoView;

    private:
        typedef DenseVector                                 DV;

    public:
        typedef flens::IndexVariable<IndexType>             IndexVariable;
        typedef densevector::ConstElementClosure<DV>        ConstElementClosure;
        typedef densevector::ElementClosure<DV>             ElementClosure;
        typedef densevector::Initializer<DV>                Initializer;

        DenseVector();

        explicit
        DenseVector(IndexType length);

        DenseVector(IndexType length, IndexType firstIndex);

        DenseVector(const Range<IndexType> &range);

        DenseVector(const Engine &engine, bool reverse=false);

        DenseVector(const DenseVector &rhs);

        template <typename RHS>
            DenseVector(const DenseVector<RHS> &rhs);

        template <typename RHS>
            DenseVector(DenseVector<RHS> &rhs);

        template <typename RHS>
            DenseVector(const Vector<RHS> &rhs);

        // -- operators --------------------------------------------------------

        Initializer
        operator=(const ElementType &value);

        DenseVector &
        operator=(const DenseVector &rhs);

        template <typename RHS>
            DenseVector &
            operator=(const Vector<RHS> &rhs);

        template <typename RHS>
            DenseVector &
            operator+=(const Vector<RHS> &rhs);

        template <typename RHS>
            DenseVector &
            operator-=(const Vector<RHS> &rhs);

        DenseVector &
        operator+=(const ElementType &rhs);

        DenseVector &
        operator-=(const ElementType &rhs);

        DenseVector &
        operator*=(const ElementType &alpha);

        DenseVector &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType index) const;

        ElementType &
        operator()(IndexType index);

        template <typename S>
            const densevector::ConstElementClosure<DenseVector,
                                                   typename Scalar<S>::Impl>
            operator()(const Scalar<S> &index) const;

        const ConstElementClosure
        operator()(const IndexVariable &index) const;

        ElementClosure
        operator()(IndexVariable &index);

        //-- views -------------------------------------------------------------

        const ConstView
        operator()(const Range<IndexType> &range) const;

        View
        operator()(const Range<IndexType> &range);

        const ConstView
        operator()(const Range<IndexType> &range,
                   IndexType firstViewIndex) const;

        View
        operator()(const Range<IndexType> &range, 
                   IndexType firstViewIndex);

        const ConstView
        operator()(const Underscore<IndexType> &all,
                   IndexType firstViewIndex) const;

        View
        operator()(const Underscore<IndexType> &all, 
                   IndexType firstViewIndex);

        const ConstView
        reverse() const;

        View
        reverse();

        // -- methods ----------------------------------------------------------
        Range<IndexType>
        range() const;

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        IndexType
        length() const;

        IndexType
        inc() const;

        IndexType
        endIndex() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        IndexType
        stride() const;

        IndexType
        indexBase() const;

        template <typename RHS>
            bool
            resize(const DenseVector<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType length,
               IndexType firstIndex = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        void
        changeIndexBase(IndexType firstIndex);

        // -- implementation ---------------------------------------------------
        const A &
        engine() const;

        A &
        engine();

        bool
        reversed() const;

    private:
        A    _array;
        bool _reverse;
};

//-- Traits --------------------------------------------------------------------

//
//  IsDenseVector
//

struct _DenseVectorChecker
{

    struct Two
    {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(DenseVector<Any>);
};

template <typename T>
struct IsDenseVector
{
    static T var;
    static const bool value = sizeof(_DenseVectorChecker::check(var))==1;
};

//
//  IsIntegerDenseVector
//

template <typename T>
struct IsIntegerDenseVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsDenseVector<TT>::value
                           && IsInteger<typename TT::ElementType>::value;
};


//
//  IsRealDenseVector
//

template <typename T>
struct IsRealDenseVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsDenseVector<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexDenseVector
//

template <typename T>
struct IsComplexDenseVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsDenseVector<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};

//
//   ImagVector
//

template <typename VZ>
struct ImagVector<DenseVector<VZ> >
{
    typedef typename ImagArray<VZ>::ConstView   ConstViewEngine;
    typedef typename ImagArray<VZ>::View        ViewEngine;
    typedef DenseVector<ConstViewEngine>        ConstView;
    typedef DenseVector<ViewEngine>             View;
};

//
//   RealVector
//

template <typename VZ>
struct RealVector<DenseVector<VZ> >
{
    typedef typename RealArray<VZ>::ConstView   ConstViewEngine;
    typedef typename RealArray<VZ>::View        ViewEngine;
    typedef DenseVector<ConstViewEngine>        ConstView;
    typedef DenseVector<ViewEngine>             View;
};

//-- DenseVector specific functions --------------------------------------------

//
//  imag
//

template <typename VZ>
    typename ImagVector<DenseVector<VZ> >::ConstView
    imag(const DenseVector<VZ> &z);

template <typename VZ>
    typename RestrictTo<IsDenseVector<VZ>::value,
             typename ImagVector<typename RemoveRef<VZ>::Type>::View>::Type
    imag(VZ &&z);


//
//  real
//

template <typename VZ>
    typename RealVector<DenseVector<VZ> >::ConstView
    real(const DenseVector<VZ> &z);

template <typename VZ>
    typename RestrictTo<IsDenseVector<VZ>::value,
             typename RealVector<typename RemoveRef<VZ>::Type>::View>::Type
    real(VZ &&z);


//
//  fillRandom
//

template <typename VX>
    typename RestrictTo<IsDenseVector<VX>::value,
             bool>::Type
    fillRandom(VX &&x);

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_DENSEVECTOR_H
