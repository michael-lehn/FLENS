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

#ifndef FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_H
#define FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/typedefs.h>
#include <flens/matrixtypes/hermitian/hermitianmatrix.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename PS>
    class SpMatrix;

template <typename PS>
    class TpMatrix;

template <typename PS>
class HpMatrix
    : public HermitianMatrix<HpMatrix<PS> >
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

        typedef HpMatrix<EngineConstView>           ConstView;
        typedef HpMatrix<EngineView>                View;
        typedef HpMatrix<EngineNoView>              NoView;

        typedef SpMatrix<EngineConstView>           ConstSymmetricView;
        typedef SpMatrix<EngineView>                SymmetricView;
        typedef SpMatrix<EngineNoView>              SymmetricNoView;

        typedef TpMatrix<EngineConstView>           ConstTriangularView;
        typedef TpMatrix<EngineView>                TriangularView;
        typedef TpMatrix<EngineNoView>              TriangularNoView;

        //-- Constructors ------------------------------------------------------

        HpMatrix(IndexType dim, StorageUpLo upLo = Upper);

        HpMatrix(const Engine &engine, StorageUpLo upLo);

        HpMatrix(const HpMatrix &rhs);

        template <typename RHS>
            HpMatrix(const HpMatrix<RHS> &rhs);

        template <typename RHS>
            HpMatrix(HpMatrix<RHS> &rhs);

        template <typename RHS,
                  class = typename RestrictTo<!IsSame<PS,RHS>::value, void>::Type>
            HpMatrix(HpMatrix<RHS> &&rhs);

        template <typename RHS>
            HpMatrix(const Matrix<RHS> &rhs);

        //-- Operators ---------------------------------------------------------

        HpMatrix &
        operator=(const ElementType &alpha);

        HpMatrix &
        operator=(const HpMatrix &rhs);

        template <typename RHS>
            HpMatrix &
            operator=(const Matrix<RHS> &rhs);

        template <typename RHS>
            HpMatrix &
            operator+=(const Matrix<RHS> &rhs);

        template <typename RHS>
            HpMatrix &
            operator-=(const Matrix<RHS> &rhs);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        HpMatrix &
        operator+=(const ElementType &alpha);

        HpMatrix &
        operator-=(const ElementType &alpha);

        HpMatrix &
        operator*=(const ElementType &alpha);

        HpMatrix &
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
            resize(const HpMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType dim,
               IndexType indexBase = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        //-- Views -------------------------------------------------------------

        // hermitian views
        const ConstView
        hermitian() const;

        View
        hermitian();

        // symmetric view
        const ConstSymmetricView
        symmetric() const;

        SymmetricView
        symmetric();

        // triangular view
        const ConstTriangularView
        triangular() const;

        TriangularView
        triangular();

        //-- Implementation ----------------------------------------------------

        const Engine &
        engine() const;

        Engine &
        engine();

    private:
        Engine       engine_;
        StorageUpLo  upLo_;
};

//-- Traits --------------------------------------------------------------------
//
//  IsHpMatrix
//
struct HpMatrixChecker_
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(AnyConversion_);

    template <typename Any>
        static char
        check(HpMatrix<Any>);
};

template <typename T>
struct IsHpMatrix
{
    static T var;
    static const bool value = sizeof(HpMatrixChecker_::check(var))==1;
};

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_H
