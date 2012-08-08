/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_H
#define FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/auxiliary/underscore.h>
#include <flens/matrixtypes/symmetric/symmetricmatrix.h>
#include <flens/typedefs.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename FS>
    class GeMatrix;

template <typename FS>
    class HeMatrix;

template <typename FS>
    class TrMatrix;

template <typename FS>
class SyMatrix
    : public SymmetricMatrix<SyMatrix<FS> >
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

        typedef GeMatrix<EngineConstView>           ConstGeneralView;
        typedef GeMatrix<EngineView>                GeneralView;
        typedef GeMatrix<EngineNoView>              GeneralNoView;

        typedef HeMatrix<EngineConstView>           ConstHermitianView;
        typedef HeMatrix<EngineView>                HermitianView;
        typedef HeMatrix<EngineNoView>              HermitianNoView;

        typedef SyMatrix<EngineConstView>           ConstView;
        typedef SyMatrix<EngineView>                View;
        typedef SyMatrix<EngineNoView>              NoView;

        typedef TrMatrix<EngineConstView>           ConstTriangularView;
        typedef TrMatrix<EngineView>                TriangularView;
        typedef TrMatrix<EngineNoView>              TriangularNoView;

        SyMatrix();

        explicit
        SyMatrix(IndexType dim);

        SyMatrix(IndexType dim, IndexType firstRow, IndexType firstCol);

        SyMatrix(const Engine &engine, StorageUpLo upLo);

        SyMatrix(const SyMatrix &rhs);

        template <typename RHS>
            SyMatrix(const SyMatrix<RHS> &rhs);

        template <typename RHS>
            SyMatrix(SyMatrix<RHS> &rhs);

        template <typename RHS>
            SyMatrix(const Matrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        void
        operator=(const ElementType &value);

        SyMatrix &
        operator=(const SyMatrix &rhs);

        template <typename RHS>
            SyMatrix &
            operator=(const Matrix<RHS> &rhs);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        // rectangular views
        const ConstGeneralView
        operator()(const Range<IndexType> &rows,
                   const Range<IndexType> &cols) const;

        GeneralView
        operator()(const Range<IndexType> &rows,
                   const Range<IndexType> &cols);

        // rectangular views (all rows selected)
        const ConstGeneralView
        operator()(const Underscore<IndexType> &,
                   const Range<IndexType> &cols) const;

        GeneralView
        operator()(const Underscore<IndexType> &,
                   const Range<IndexType> &cols);

        // rectangular views (all columns selected)
        const ConstGeneralView
        operator()(const Range<IndexType> &rows,
                   const Underscore<IndexType> &) const;

        GeneralView
        operator()(const Range<IndexType> &rows,
                   const Underscore<IndexType> &);

        // row view (vector view)
        const ConstVectorView
        operator()(IndexType row, const Underscore<IndexType> &) const;

        VectorView
        operator()(IndexType row, const Underscore<IndexType> &);

        const ConstVectorView
        operator()(IndexType row, const Range<IndexType> &cols) const;

        VectorView
        operator()(IndexType row, const Range<IndexType> &cols);

        // column view (vector view)
        const ConstVectorView
        operator()(const Underscore<IndexType> &, IndexType col) const;

        VectorView
        operator()(const Underscore<IndexType> &, IndexType col);

        const ConstVectorView
        operator()(const Range<IndexType> &rows, IndexType col) const;

        VectorView
        operator()(const Range<IndexType> &rows, IndexType col);

        // -- views ------------------------------------------------------------

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

        // triangular views
        const ConstTriangularView
        triangular() const;

        TriangularView
        triangular();

        // diag views
        const ConstVectorView
        diag(IndexType d) const;

        VectorView
        diag(IndexType d);

        // -- methods ----------------------------------------------------------

        IndexType
        dim() const;

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

        StorageOrder
        order() const;

        template <typename RHS>
            bool
            resize(const SyMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType dim,
               IndexType firstIndex = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        // -- implementation ---------------------------------------------------

        const Engine &
        engine() const;

        Engine &
        engine();

        StorageUpLo
        upLo() const;

        StorageUpLo &
        upLo();

    private:
        Engine      _engine;
        StorageUpLo _upLo;
};

//-- Traits --------------------------------------------------------------------
//
//  IsSyMatrix
//
struct _SyMatrixChecker
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(_AnyConversion);

    template <typename Any>
        static char
        check(SyMatrix<Any>);
};

template <typename T>
struct IsSyMatrix
{
    static T var;
    static const bool value = sizeof(_SyMatrixChecker::check(var))==1;
};

//
//  IsRealSyMatrix
//
template <typename T>
struct IsRealSyMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSyMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexSyMatrix
//
template <typename T>
struct IsComplexSyMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsSyMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};

} // namespace flens

#endif // FLENS_MATRIXTYPES_SYMMETRIC_IMPL_SYMATRIX_H
