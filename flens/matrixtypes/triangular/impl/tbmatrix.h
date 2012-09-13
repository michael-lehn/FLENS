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

#ifndef FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_H
#define FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/triangular/triangularmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename FS>
    class GbMatrix;

template <typename FS>
    class HbMatrix;

template <typename FS>
    class SbMatrix;

template <typename FS>
class TbMatrix
    : public TriangularMatrix<TbMatrix<FS> >
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
        typedef DenseVector<ConstArrayView>         ConstVectorView;
        typedef DenseVector<ArrayView>              VectorView;
        typedef DenseVector<Array>                  Vector;

        typedef GbMatrix<EngineConstView>           ConstGeneralView;
        typedef GbMatrix<EngineView>                GeneralView;
        typedef GbMatrix<EngineNoView>              GeneralNoView;

        typedef HbMatrix<EngineConstView>           ConstHermitianView;
        typedef HbMatrix<EngineView>                HermitianView;
        typedef HbMatrix<EngineNoView>              HermitianNoView;

        typedef SbMatrix<EngineConstView>           ConstSymmetricView;
        typedef SbMatrix<EngineView>                SymmetricView;
        typedef SbMatrix<EngineNoView>              SymmetricNoView;

        typedef TbMatrix<EngineConstView>           ConstView;
        typedef TbMatrix<EngineView>                View;
        typedef TbMatrix<EngineNoView>              NoView;

        TbMatrix();

        explicit
        TbMatrix(IndexType dim, StorageUpLo upLo, IndexType numOffDiags,
                 Diag diag = NonUnit);

        TbMatrix(const Engine &engine, StorageUpLo upLo, Diag diag = NonUnit);

        TbMatrix(const TbMatrix &rhs);

        template <typename RHS>
            TbMatrix(const TbMatrix<RHS> &rhs);

        template <typename RHS>
            TbMatrix(TbMatrix<RHS> &rhs);

        template <typename RHS>
            TbMatrix(const Matrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        TbMatrix &
        operator=(const TbMatrix &rhs);

        template <typename RHS>
            TbMatrix &
            operator=(const Matrix<RHS> &rhs);

        template <typename RHS>
            TbMatrix &
            operator+=(const Matrix<RHS> &rhs);

        template <typename RHS>
            TbMatrix &
            operator-=(const Matrix<RHS> &rhs);

        TbMatrix<FS> &
        operator=(const ElementType &alpha);

        TbMatrix<FS> &
        operator+=(const ElementType &alpha);

        TbMatrix<FS> &
        operator-=(const ElementType &alpha);

        TbMatrix &
        operator*=(const ElementType &alpha);

        TbMatrix &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // general views
        const ConstGeneralView
        general() const;

        GeneralView
        general();

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

        // triangular view
        const ConstView
        triangular() const;

        View
        triangular();

        VectorView
        viewDiag(IndexType diag = IndexType(0));

        const ConstVectorView
        viewDiag(IndexType diag = IndexType(0)) const;

        // -- views ------------------------------------------------------------


        // -- methods ----------------------------------------------------------
        IndexType
        dim() const;

        IndexType
        numCols() const;

        IndexType
        numRows() const;

        IndexType
        numOffDiags() const;

        IndexType
        numSubDiags() const;

        IndexType
        numSuperDiags() const;

        // for element access
        IndexType
        firstIndex() const;

        IndexType
        lastIndex() const;

        IndexType
        leadingDimension() const;

        StorageOrder
        order() const;

        const ElementType *
        data() const;

        ElementType *
        data();

        template <typename RHS>
            bool
            resize(const TbMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType dim,
               IndexType numOffDiags,
               IndexType firstColIndex = Engine::defaultIndexBase,
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

        StorageUpLo &
        upLo();

        Diag
        diag() const;

        Diag &
        diag();

    private:
        Engine       _engine;
        StorageUpLo  _upLo;
        Diag         _diag;
};

//-- Traits --------------------------------------------------------------------
//
//  IsTbMatrix
//
struct _TbMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(TbMatrix<Any>);
};

template <typename T>
struct IsTbMatrix
{
    static T var;
    static const bool value = sizeof(_TbMatrixChecker::check(var))==1;
};

//
//  IsRealTbMatrix
//
template <typename T>
struct IsRealTbMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTbMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexTbMatrix
//
template <typename T>
struct IsComplexTbMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsTbMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};


} // namespace flens

#endif // FLENS_MATRIXTYPES_TRIANGULAR_IMPL_TBMATRIX_H
