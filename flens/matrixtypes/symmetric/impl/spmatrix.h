/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_H
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/symmetric/symmetricmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename FS>
    class HpMatrix;

template <typename FS>
    class TpMatrix;

template <typename FS>
class SpMatrix
    : public SymmetricMatrix<SpMatrix<FS> >
{
    public:
        typedef FS                                  Engine;
        typedef typename Engine::ElementType        ElementType;
        typedef typename Engine::IndexType          IndexType;

        // view types from Engine
        typedef typename Engine::ConstView          EngineConstView;
        typedef typename Engine::View               EngineView;
        typedef typename Engine::NoView             EngineNoView;

        typedef typename Engine::ConstArrayView     ConstArrayView;
        typedef typename Engine::ArrayView          ArrayView;
        typedef typename Engine::Array              Array;

        // view types for
        // view types for
        typedef DenseVector<ConstArrayView>         ConstVectorView;
        typedef DenseVector<ArrayView>              VectorView;
        typedef DenseVector<Array>                  Vector;

        typedef HpMatrix<EngineConstView>           ConstHermitianView;
        typedef HpMatrix<EngineView>                HermitianView;
        typedef HpMatrix<EngineNoView>              HermitianNoView;

        typedef SpMatrix<EngineConstView>           ConstView;
        typedef SpMatrix<EngineView>                View;
        typedef SpMatrix<EngineNoView>              NoView;

        typedef TpMatrix<EngineConstView>           ConstTriangularView;
        typedef TpMatrix<EngineView>                TriangularView;
        typedef TpMatrix<EngineNoView>              TriangularNoView;

        SpMatrix();

        explicit
        SpMatrix(IndexType dim);

        SpMatrix(const Engine &engine);

        SpMatrix(const SpMatrix &rhs);

        template <typename RHS>
            SpMatrix(const SpMatrix<RHS> &rhs);

        template <typename RHS>
            SpMatrix(SpMatrix<RHS> &rhs);

        template <typename RHS>
            SpMatrix(const Matrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        SpMatrix &
        operator=(const SpMatrix &rhs);

        template <typename RHS>
            SpMatrix &
            operator=(const Matrix<RHS> &rhs);

        SpMatrix<FS> &
        operator=(const ElementType &alpha);

        SpMatrix<FS> &
        operator+=(const ElementType &alpha);

        SpMatrix<FS> &
        operator-=(const ElementType &alpha);

        SpMatrix<FS> &
        operator*=(const ElementType &alpha);

        SpMatrix<FS> &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // -- views ------------------------------------------------------------

        // hermitian views
        const ConstHermitianView
        hermitian() const;

        HermitianView
        hermitian();

        // symmetric views
        const ConstView
        symmetric() const;

        View
        symmetric();

        // triangular views
        const ConstTriangularView
        triangular() const;

        TriangularView
        triangular();

        // -- methods ----------------------------------------------------------

        IndexType
        dim() const;

        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        StorageOrder
        order() const;

        template <typename RHS>
            bool
            resize(const SpMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType dim,
               IndexType firstIndex = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        bool
        fillRandom();

        // -- implementation ---------------------------------------------------

        const Engine &
        engine() const;

        Engine &
        engine();

        StorageUpLo
        upLo() const;

    private:
        Engine      _engine;
        StorageUpLo _upLo;
};

//-- Traits --------------------------------------------------------------------
//
//  IsSpMatrix
//
struct _SpMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(SpMatrix<Any>);
};

template <typename T>
struct IsSpMatrix
{
    static T var;
    static const bool value = sizeof(_SpMatrixChecker::check(var))==1;
};

//
//  IsRealSpMatrix
//
template <typename T>
struct IsRealSpMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSpMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexSpMatrix
//
template <typename T>
struct IsComplexSpMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSpMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};


} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SPMATRIX_H
