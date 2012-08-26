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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCOORDMATRIX_H
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCOORDMATRIX_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/symmetric/symmetricmatrix.h>
#include <flens/typedefs.h>

namespace flens {

template <typename CS>
class SyCoordMatrix
    : public SymmetricMatrix<SyCoordMatrix<CS> >
{
    public:
        typedef CS                              Engine;
        typedef typename Engine::ElementType    ElementType;
        typedef typename Engine::IndexType      IndexType;
        typedef typename Engine::ElementProxy   ElementProxy;

        // -- constructors -----------------------------------------------------
        explicit
        SyCoordMatrix(IndexType   dim,
                      StorageUpLo upLo = Upper,
                      IndexType   densityEstimate = 1,
                      IndexType   indexBase = CS::defaultIndexBase);

        // -- operators --------------------------------------------------------
        ElementProxy
        operator()(IndexType row, IndexType col);

        // -- methods ----------------------------------------------------------
        IndexType
        dim() const;

        IndexType
        numRows() const;

        IndexType
        numCols() const;

        IndexType
        indexBase() const;

        IndexType
        firstRow() const;

        IndexType
        lastRow() const;

        IndexType
        firstCol() const;

        IndexType
        lastCol() const;

        StorageUpLo
        upLo() const;

        // -- implementation ---------------------------------------------------
        const Engine &
        engine() const;

    private:
        // forbidden:
        SyCoordMatrix(const SyCoordMatrix &rhs);

        Engine       _engine;
        StorageUpLo  _upLo;
};

//-- Traits --------------------------------------------------------------------
//
//  IsSyCoordMatrix
//
struct _SyCoordMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(const SyCoordMatrix<Any> &);
};

template <typename T>
struct IsSyCoordMatrix
{
    static T var;
    static const bool value = sizeof(_SyCoordMatrixChecker::check(var))==1;
};

//
//  IsRealSyCoordMatrix
//
template <typename T>
struct IsRealCoordSyMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSyCoordMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexSyCoordMatrix
//
template <typename T>
struct IsComplexSyCoordMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSyCoordMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};


} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYCOORDMATRIX_H

