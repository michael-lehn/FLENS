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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_H
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/general/generalmatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename CS>
class GeCoordMatrix
    : public GeneralMatrix<GeCoordMatrix<CS> >
{
    public:
        typedef CS                              Engine;
        typedef typename Engine::ElementType    ElementType;
        typedef typename Engine::IndexType      IndexType;
        typedef typename Engine::ElementProxy   ElementProxy;

        // -- constructors -----------------------------------------------------
        GeCoordMatrix(IndexType numRows, IndexType numCols,
                      IndexType densityEstimate = 1);

        GeCoordMatrix(IndexType numRows, IndexType numCols,
                      IndexType firstRow, IndexType firstCol,
                      IndexType densityEstimate = 1);

        // -- operators --------------------------------------------------------
        ElementProxy
        operator()(IndexType row, IndexType col);

        // -- methods ----------------------------------------------------------
        const IndexType
        numRows() const;

        const IndexType
        numCols() const;

        const IndexType
        firstRow() const;

        const IndexType
        lastRow() const;

        const IndexType
        firstCol() const;

        const IndexType
        lastCol() const;

        // -- implementation ---------------------------------------------------
        const Engine &
        engine() const;

    private:
        // forbidden:
        GeCoordMatrix(const GeCoordMatrix &rhs);

        Engine _engine;
};

//-- Traits --------------------------------------------------------------------
//
//  IsGeCoordMatrix
//
struct _GeCoordMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(const GeCoordMatrix<Any> &);
};

template <typename T>
struct IsGeCoordMatrix
{
    static T var;
    static const bool value = sizeof(_GeCoordMatrixChecker::check(var))==1;
};

//
//  IsRealGeCoordMatrix
//
template <typename T>
struct IsRealCoordGeMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeCoordMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexGeCoordMatrix
//
template <typename T>
struct IsComplexGeCoordMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeCoordMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};


} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GECOORDMATRIX_H

