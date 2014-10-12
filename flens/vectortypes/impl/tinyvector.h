/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef FLENS_VECTORTYPES_IMPL_TINYVECTOR_H
#define FLENS_VECTORTYPES_IMPL_TINYVECTOR_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/vectortypes/vector.h>

namespace flens {

template <typename TA>
class TinyVector
    : public Vector<TinyVector<TA> >
{
    public:
        typedef TA                                          Engine;
        typedef typename Engine::ElementType                ElementType;
        typedef typename Engine::IndexType                  IndexType;

        // -- constructors -----------------------------------------------------
        TinyVector();

        TinyVector(const Engine &engine);

        template <typename RHS>
            TinyVector(const TinyVector<RHS> &rhs);

        template <typename RHS>
            TinyVector(TinyVector<RHS> &rhs);

        template <typename RHS>
            TinyVector(const Vector<RHS> &rhs);

        // -- operators --------------------------------------------------------

        TinyVector &
        operator=(const TinyVector &rhs);

        template <typename RHS>
            TinyVector &
            operator=(const Vector<RHS> &rhs);

        template <typename RHS>
            TinyVector &
            operator+=(const Vector<RHS> &rhs);

        template <typename RHS>
            TinyVector &
            operator-=(const Vector<RHS> &rhs);

        TinyVector &
        operator+=(const ElementType &rhs);

        TinyVector &
        operator-=(const ElementType &rhs);

        TinyVector &
        operator*=(const ElementType &alpha);

        TinyVector &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType index) const;

        ElementType &
        operator()(IndexType index);

        // -- methods ----------------------------------------------------------

        IndexType
        length() const;

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        IndexType
        stride() const;

        void
        fill(const ElementType &value = ElementType(0));

        // -- implementation ---------------------------------------------------
        const TA &
        engine() const;

        TA &
        engine();

    private:
        TA    array_;
};

//-- Traits --------------------------------------------------------------------

//
//  IsTinyVector
//

struct TinyVectorChecker_
{

    struct Two
    {
        char x;
        char y;
    };

    static Two
    check(AnyConversion_);

    template <typename Any>
        static char
        check(TinyVector<Any>);
};

template <typename T>
struct IsTinyVector
{
    static T var;
    static const bool value = sizeof(TinyVectorChecker_::check(var))==1;
};

//
//  IsIntegerTinyVector
//

template <typename T>
struct IsIntegerTinyVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTinyVector<TT>::value
                           && IsInteger<typename TT::ElementType>::value;
};


//
//  IsRealTinyVector
//

template <typename T>
struct IsRealTinyVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTinyVector<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexTinyVector
//

template <typename T>
struct IsComplexTinyVector
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTinyVector<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};

} // namespace flens

#endif // FLENS_VECTORTYPES_IMPL_TINYVECTOR_H
