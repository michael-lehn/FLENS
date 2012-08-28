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
#include <flens/typedefs.h>
#include <flens/matrixtypes/hermitian/hermitianmatrix.h>

namespace flens {

// forward declarations
template <typename FS>
    class SpMatrix;

template <typename FS>
    class TpMatrix;

template <typename FS>
class HpMatrix
    : public HermitianMatrix<HpMatrix<FS> >
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

        HpMatrix();

        explicit
        HpMatrix(IndexType dim);

        HpMatrix(IndexType dim,
                 IndexType firstIndex);

        HpMatrix(const Engine &engine);

        HpMatrix(const HpMatrix &rhs);

        template <typename RHS>
            HpMatrix(const HpMatrix<RHS> &rhs);

        template <typename RHS>
            HpMatrix(HpMatrix<RHS> &rhs);

        template <typename RHS>
            HpMatrix(const Matrix<RHS> &rhs);


        // -- operators --------------------------------------------------------
        HpMatrix &
        operator=(const HpMatrix &rhs);

        template <typename RHS>
            HpMatrix &
            operator=(const Matrix<RHS> &rhs);

        HpMatrix<FS> &
        operator=(const ElementType &alpha);

        HpMatrix<FS> &
        operator+=(const ElementType &alpha);

        HpMatrix<FS> &
        operator-=(const ElementType &alpha);

        HpMatrix<FS> &
        operator*=(const ElementType &alpha);

        HpMatrix<FS> &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // -- views ------------------------------------------------------------

        // general views
        ConstView
        hermitian() const;

        View
        hermitian();

        // symmetric view
        ConstSymmetricView
        symmetric() const;

        SymmetricView
        symmetric();

        // triangular view
        ConstTriangularView
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
            resize(const HpMatrix<RHS> &rhs,
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
};

//-- Traits --------------------------------------------------------------------
//
//  IsHpMatrix
//
struct _HpMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(HpMatrix<Any>);
};

template <typename T>
struct IsHpMatrix
{
    static T var;
    static const bool value = sizeof(_HpMatrixChecker::check(var))==1;
};

} // namespace flens

#endif // FLENS_MATRIXTYPES_HERMITIAN_IMPL_HPMATRIX_H
