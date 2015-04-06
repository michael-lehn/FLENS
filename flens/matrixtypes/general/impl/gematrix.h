/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_H
#define FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_H 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/general/generalmatrix.h>
#include <flens/matrixtypes/general/impl/ge/constelementclosure.h>
#include <flens/matrixtypes/general/impl/ge/elementclosure.h>
#include <flens/matrixtypes/general/impl/ge/initializer.h>
#include <flens/matrixtypes/general/impl/imagmatrixclosure.h>
#include <flens/matrixtypes/general/impl/realmatrixclosure.h>
#include <flens/scalartypes/impl/indexvariable.h>
#include <flens/typedefs.h>

namespace flens {

// forward declarations
template <typename A>
    class DenseVector;

template <typename FS>
    class HeMatrix;

template <typename FS>
    class SyMatrix;

template <typename FS>
    class TrMatrix;

template <typename FS>
class GeMatrix
    : public GeneralMatrix<GeMatrix<FS> >
{
    public:
        typedef FS                                      Engine;
        typedef typename Engine::ElementType            ElementType;
        typedef typename Engine::IndexType              IndexType;

        // view types from Engine
        typedef typename Engine::ConstView              EngineConstView;
        typedef typename Engine::View                   EngineView;
        typedef typename Engine::NoView                 EngineNoView;

        typedef typename Engine::ConstArrayView         ConstArrayView;
        typedef typename Engine::ArrayView              ArrayView;
        typedef typename Engine::Array                  Array;

        // view types
        typedef DenseVector<ConstArrayView>             ConstVectorView;
        typedef DenseVector<ArrayView>                  VectorView;
        typedef DenseVector<Array>                      Vector;

        typedef GeMatrix<EngineConstView>               ConstView;
        typedef GeMatrix<EngineView>                    View;
        typedef GeMatrix<EngineNoView>                  NoView;

        typedef HeMatrix<EngineConstView>               ConstHermitianView;
        typedef HeMatrix<EngineView>                    HermitianView;
        typedef HeMatrix<EngineNoView>                  HermitianNoView;

        typedef SyMatrix<EngineConstView>               ConstSymmetricView;
        typedef SyMatrix<EngineView>                    SymmetricView;
        typedef SyMatrix<EngineNoView>                  SymmetricNoView;

        typedef TrMatrix<EngineConstView>               ConstTriangularView;
        typedef TrMatrix<EngineView>                    TriangularView;
        typedef TrMatrix<EngineNoView>                  TriangularNoView;

    private:
        typedef GeMatrix                                GE;

    public:
        typedef flens::IndexVariable<IndexType>         IndexVariable;
        typedef gematrix::ConstElementClosure<GE>       ConstElementClosure;
        typedef gematrix::ElementClosure<GE>            ElementClosure;
        typedef gematrix::Initializer<GE>               Initializer;

        // -- constructors -----------------------------------------------------
        GeMatrix();

        GeMatrix(IndexType numRows, IndexType numCols);

        GeMatrix(IndexType numRows, IndexType numCols,
                 IndexType firstRow, IndexType firstCol);

        GeMatrix(const Range<IndexType> &rowRange,
                 const Range<IndexType> &colRange);

        GeMatrix(const Engine &engine);

        GeMatrix(const GeMatrix &rhs);

        template <typename RHS>
            GeMatrix(const GeMatrix<RHS> &rhs);

        template <typename RHS>
            GeMatrix(GeMatrix<RHS> &rhs);

        template <typename RHS,
                  class = typename RestrictTo<!IsSame<FS,RHS>::value, void>::Type>
            GeMatrix(GeMatrix<RHS> &&rhs);

        template <typename RHS>
            GeMatrix(const Matrix<RHS> &rhs);

        template <typename VECTOR>
            GeMatrix(IndexType numRows, IndexType numCols,
                     VECTOR &&rhs);

        template <typename VECTOR>
            GeMatrix(IndexType numRows, IndexType numCols,
                     VECTOR &&rhs,
                     IndexType leadingDimension);

        // -- operators --------------------------------------------------------
        Initializer
        operator=(const ElementType &value);

        GeMatrix &
        operator=(const GeMatrix &rhs);

        template <typename RHS>
            GeMatrix &
            operator=(const Matrix<RHS> &rhs);

        template <typename RHS>
            GeMatrix &
            operator+=(const Matrix<RHS> &rhs);

        template <typename RHS>
            GeMatrix &
            operator-=(const Matrix<RHS> &rhs);

        GeMatrix &
        operator+=(const ElementType &alpha);

        GeMatrix &
        operator-=(const ElementType &alpha);

        GeMatrix &
        operator*=(const ElementType &alpha);

        GeMatrix &
        operator/=(const ElementType &alpha);

        const ElementType &
        operator()(IndexType row, IndexType col) const;

        ElementType &
        operator()(IndexType row, IndexType col);

        template <typename S>
            const gematrix::ConstElementClosure<GeMatrix,
                                                typename Scalar<S>::Impl>
            operator()(const Scalar<S> &row, const Scalar<S> &col) const;

        const ConstElementClosure
        operator()(const IndexVariable &row, const IndexVariable &col) const;

        ElementClosure
        operator()(IndexVariable &row, IndexVariable &col);

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

        Range<IndexType>
        rows() const;

        Range<IndexType>
        cols() const;

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
            resize(const GeMatrix<RHS> &rhs,
                   const ElementType &value = ElementType());

        bool
        resize(IndexType numRows, IndexType numCols,
               IndexType firstRowIndex = Engine::defaultIndexBase,
               IndexType firstColIndex = Engine::defaultIndexBase,
               const ElementType &value = ElementType());

        bool
        fill(const ElementType &value = ElementType(0));

        void
        changeIndexBase(IndexType firstRowIndex, IndexType firstColIndex);

        // -- views ------------------------------------------------------------
        // vectorize matrix
        const ConstVectorView
        vectorView() const;

        VectorView
        vectorView();

        // vectorize matrix and select range
        const ConstVectorView
        vectorView(IndexType from, IndexType to) const;

        VectorView
        vectorView(IndexType from, IndexType to);

        // diag views
        const ConstVectorView
        diag(IndexType d) const;

        VectorView
        diag(IndexType d);

        // anti diag views
        const ConstVectorView
        antiDiag(IndexType d) const;

        VectorView
        antiDiag(IndexType d);

        // triangular views
        const ConstTriangularView
        upper() const;

        TriangularView
        upper();

        const ConstTriangularView
        upperUnit() const;

        TriangularView
        upperUnit();

        const ConstTriangularView
        strictUpper() const;

        TriangularView
        strictUpper();

        const ConstTriangularView
        lower() const;

        TriangularView
        lower();

        const ConstTriangularView
        lowerUnit() const;

        TriangularView
        lowerUnit();

        const ConstTriangularView
        strictLower() const;

        TriangularView
        strictLower();

        // rectangular views
        const ConstView
        operator()(const Range<IndexType> &rows,
                   const Range<IndexType> &cols) const;

        View
        operator()(const Range<IndexType> &rows,
                   const Range<IndexType> &cols);

        template <typename RHS>
            const ConstView
            operator()(const GeMatrix<RHS> &A) const;

        template <typename RHS>
            View
            operator()(const GeMatrix<RHS> &A);

        // rectangular views (all rows selected)
        const ConstView
        operator()(const Underscore<IndexType> &,
                   const Range<IndexType> &cols) const;

        View
        operator()(const Underscore<IndexType> &,
                   const Range<IndexType> &cols);

        // rectangular views (all columns selected)
        const ConstView
        operator()(const Range<IndexType> &rows,
                   const Underscore<IndexType> &) const;

        View
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
struct GeMatrixChecker_
{

    struct Two {
        char x;
        char y;
    };

    static Two
    check(AnyConversion_);

    template <typename Any>
        static char
        check(GeMatrix<Any>);
};

template <typename T>
struct IsGeMatrix
{
    static T var;
    static const bool value = sizeof(GeMatrixChecker_::check(var))==1;
};

//
//  IsRealGeMatrix
//
template <typename T>
struct IsRealGeMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeMatrix<TT>::value
                           && IsNotComplex<typename TT::ElementType>::value;
};

//
//  IsComplexGeMatrix
//
template <typename T>
struct IsComplexGeMatrix
{
    typedef typename std::remove_reference<T>::type  TT;

    static const bool value = IsGeMatrix<TT>::value
                           && IsComplex<typename TT::ElementType>::value;
};

//-- GeMatrix specific functions -----------------------------------------------

//
//  imag
//
template <typename MZ>
    ImagConstMatrixClosure<GeMatrix<MZ> >
    imag(const GeMatrix<MZ> &Z);

template <typename MZ>
    typename RestrictTo<IsGeMatrix<MZ>::value,
             ImagMatrixClosure<MZ> >::Type
    imag(MZ &&Z);

//
//  real
//
template <typename MZ>
    RealConstMatrixClosure<GeMatrix<MZ> >
    real(const GeMatrix<MZ> &Z);

template <typename MZ>
    typename RestrictTo<IsGeMatrix<MZ>::value,
             RealMatrixClosure<MZ> >::Type
    real(MZ &&Z);

//
//  fillRandom
//

template <typename MA>
    typename RestrictTo<IsGeMatrix<MA>::value,
             bool>::Type
    fillRandom(MA &&A);

} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_GEMATRIX_H
