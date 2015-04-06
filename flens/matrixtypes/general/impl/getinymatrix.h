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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_H
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/general/generalmatrix.h>

namespace flens {

template <typename TFS>
class GeTinyMatrix
    : public GeneralMatrix<GeTinyMatrix<TFS> >
{
    public:
        typedef TFS                                     Engine;
        typedef typename Engine::ElementType            ElementType;
        typedef typename Engine::IndexType              IndexType;

        // -- constructors -----------------------------------------------------
        GeTinyMatrix();

        GeTinyMatrix(const Engine &engine);

       template <typename RHS>
            GeTinyMatrix(const GeTinyMatrix<RHS> &rhs);

        template <typename RHS>
            GeTinyMatrix(GeTinyMatrix<RHS> &rhs);

        template <typename RHS,
                  class = typename RestrictTo<!IsSame<TFS,RHS>::value, void>::Type>
            GeTinyMatrix(GeTinyMatrix<RHS> &&rhs);

        template <typename RHS>
            GeTinyMatrix(const Matrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        GeTinyMatrix &
        operator=(const GeTinyMatrix &rhs);

        template <typename RHS>
            GeTinyMatrix &
            operator=(const Matrix<RHS> &rhs);

        template <typename RHS>
            GeTinyMatrix &
            operator+=(const Matrix<RHS> &rhs);

        template <typename RHS>
            GeTinyMatrix &
            operator-=(const Matrix<RHS> &rhs);

        GeTinyMatrix &
        operator+=(const ElementType &alpha);

        GeTinyMatrix &
        operator-=(const ElementType &alpha);

        GeTinyMatrix &
        operator*=(const ElementType &alpha);

        GeTinyMatrix &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // -- methods ----------------------------------------------------------
        IndexType
        numRows() const;

        IndexType
        numCols() const;

        IndexType
        firstRow() const;

        IndexType
        lastRow() const;

        IndexType
        firstCol() const;

        IndexType
        lastCol() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        IndexType
        leadingDimension() const;

        void
        fill(const ElementType &value = ElementType(0));

        // -- implementation ---------------------------------------------------
        const Engine &
        engine() const;

        Engine &
        engine();

    private:
        Engine engine_;
};

//-- Traits --------------------------------------------------------------------
//
//  IsGeMatrix
//
struct GeTinyMatrixChecker_
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(AnyConversion_);

    template <typename Any>
        static char
        check(GeTinyMatrix<Any>);
};

template <typename T>
struct IsGeTinyMatrix
{
    static T var;
    static const bool value = sizeof(GeTinyMatrixChecker_::check(var))==1;
};

//
//  IsRealGeTinyMatrix
//
template <typename T>
struct IsRealGeTinyMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeTinyMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexGeTinyMatrix
//
template <typename T>
struct IsComplexGeTinyMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeTinyMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GETINYMATRIX_H
