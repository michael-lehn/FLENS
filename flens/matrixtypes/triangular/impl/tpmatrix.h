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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_H
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/triangular/triangularmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename PS>
    class HpMatrix;

template <typename PS>
    class SpMatrix;

template <typename PS>
class TpMatrix
    : public TriangularMatrix<TpMatrix<PS> >
{
    public:
        typedef PS                                  Engine;
        typedef typename Engine::ElementType        ElementType;
        typedef typename Engine::IndexType          IndexType;

        // view types from Engine
        typedef typename Engine::ConstView          EngineConstView;
        typedef typename Engine::View               EngineView;
        typedef typename Engine::NoView             EngineNoView;

        typedef typename Engine::ConstArrayView     ConstArrayView;
        typedef typename Engine::ArrayView          ArrayView;
        typedef typename Engine::Array              Array;

        // view types
        typedef DenseVector<ConstArrayView>         ConstVectorView;
        typedef DenseVector<ArrayView>              VectorView;
        typedef DenseVector<Array>                  Vector;

        typedef HpMatrix<EngineConstView>           ConstHermitianView;
        typedef HpMatrix<EngineView>                HermitianView;
        typedef HpMatrix<EngineNoView>              HermitianNoView;

        typedef SpMatrix<EngineConstView>           ConstSymmetricView;
        typedef SpMatrix<EngineView>                SymmetricView;
        typedef SpMatrix<EngineNoView>              SymmetricNoView;

        typedef TpMatrix<EngineConstView>           ConstView;
        typedef TpMatrix<EngineView>                View;
        typedef TpMatrix<EngineNoView>              NoView;

        //-- Constructors ------------------------------------------------------

        TpMatrix(IndexType dim, StorageUpLo upLo, Diag diag = NonUnit);

        TpMatrix(const Engine &engine, StorageUpLo upLo, Diag diag = NonUnit);

        TpMatrix(const TpMatrix &rhs);

        template <typename RHS>
            TpMatrix(const TpMatrix<RHS> &rhs);

        template <typename RHS>
            TpMatrix(TpMatrix<RHS> &rhs);

        template <typename RHS>
            TpMatrix(const Matrix<RHS> &rhs);

        //-- Operators ---------------------------------------------------------

        TpMatrix &
        operator=(const ElementType &alpha);

        TpMatrix &
        operator=(const TpMatrix &rhs);

        template <typename RHS>
            TpMatrix &
            operator=(const Matrix<RHS> &rhs);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        TpMatrix &
        operator+=(const ElementType &alpha);

        TpMatrix &
        operator-=(const ElementType &alpha);

        TpMatrix &
        operator*=(const ElementType &alpha);

        TpMatrix &
        operator/=(const ElementType &alpha);

        //-- Methods -----------------------------------------------------------

        IndexType
        numRows() const;

        IndexType
        numCols() const;

        IndexType
        dim() const;

        IndexType
        firstRow() const;

        IndexType
        lastRow() const;

        IndexType
        firstCol() const;

        IndexType
        lastCol() const;

        IndexType
        indexBase() const;

        StorageUpLo
        upLo() const;

        StorageUpLo &
        upLo();

        Diag
        diag() const;

        Diag &
        diag();

        const ElementType *
        data() const;

        ElementType *
        data();

        StorageOrder
        order() const;

        bool
        fill(const ElementType &value = ElementType(0));

        template <typename RHS>
            bool
            resize(const TpMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType dim,
               IndexType indexBase = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        //-- Views -------------------------------------------------------------

        // hermitian views
        const ConstHermitianView
        hermitian() const;

        HermitianView
        hermitian();

        // symmetric views
        const ConstSymmetricView
        symmetric() const;

        SymmetricView
        symmetric();

        // triangular views
        const ConstView
        triangular() const;

        View
        triangular();

        //-- Implementation ----------------------------------------------------

        const Engine &
        engine() const;

        Engine &
        engine();

    private:
        Engine       engine_;
        StorageUpLo  upLo_;
        Diag         diag_;
};

//-- Traits --------------------------------------------------------------------
//
//  IsTpMatrix
//
struct TpMatrixChecker_
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(AnyConversion_);

    template <typename Any>
        static char
        check(TpMatrix<Any>);
};

template <typename T>
struct IsTpMatrix
{
    static T var;
    static const bool value = sizeof(TpMatrixChecker_::check(var))==1;
};

//
//  IsRealTpMatrix
//
template <typename T>
struct IsRealTpMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTpMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexTpMatrix
//
template <typename T>
struct IsComplexTpMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTpMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};


} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TPMATRIX_H
