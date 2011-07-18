/*
 *   Copyright (c) 2009, Michael Lehn
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

#ifndef FLENS_HERMITIANMATRIX_H
#define FLENS_HERMITIANMATRIX_H 1

#include <flens/array.h>
#include <flens/densevector.h>
#include <flens/matvec.h>
#include <flens/storage.h>
#include <flens/triangularmatrix.h>

namespace flens {

// == HeMatrix =================================================================

template <typename FS>
class HeMatrix
    : public HermitianMatrix<HeMatrix<FS> >
{
    public:
        typedef typename HeMatrix<FS>::ElementType T;

        // view types from FS
        typedef typename FS::ConstView          ConstFSView;
        typedef typename FS::View               FSView;
        typedef typename FS::View               FSNoView;

        // view types for HeMatrix
        typedef GeMatrix<ConstFSView>           ConstGeneralView;
        typedef GeMatrix<FSView>                GeneralView;
        typedef GeMatrix<FSNoView>              GeneralNoView;

        typedef TrMatrix<ConstFSView>           ConstTriangularView;
        typedef TrMatrix<FSView>                TriangularView;
        typedef TrMatrix<FSNoView>              TriangularNoView;

        HeMatrix();

        HeMatrix(int dim, StorageUpLo upLo, int firstIndex=1);

        HeMatrix(const FS &fs, StorageUpLo upLo);

        HeMatrix(const HeMatrix<FS> &rhs);

        template <typename RHS>
            HeMatrix(const HeMatrix<RHS> &rhs);

        HeMatrix(const TrMatrix<FS> &rhs);

        template <typename RHS>
            HeMatrix(const TrMatrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        HeMatrix &
        operator*=(T alpha);

        HeMatrix &
        operator/=(T alpha);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        // -- views ------------------------------------------------------------

        // general views
        ConstGeneralView
        general() const;

        GeneralView
        general();

        // triangular views
        ConstTriangularView
        triangular() const;

        TriangularView
        triangular();

        // -- methods ----------------------------------------------------------

        // for BLAS/LAPACK
        StorageUpLo
        upLo() const;

        int
        dim() const;

        int
        leadingDimension() const;

        const T *
        data() const;

        T *
        data();

        // for element access
        int
        firstRow() const;

        int
        lastRow() const;

        int
        firstCol() const;

        int
        lastCol() const;

        Range
        rows() const;

        Range
        cols() const;

        // -- implementation ---------------------------------------------------

        const FS &
        engine() const;

        FS &
        engine();

    private:
        FS _fs;
        StorageUpLo _upLo;
};

template <typename FS>
struct TypeInfo<HeMatrix<FS> >
{
    typedef HeMatrix<FS> Impl;
    typedef typename FS::ElementType ElementType;
};

// == HbMatrix =================================================================

template <typename BS>
    class GbMatrix;

template <typename BS>
    class TbMatrix;

template <typename BS>
class HbMatrix
    : public HermitianMatrix<HbMatrix<BS> >
{
    public:
        typedef typename HbMatrix<BS>::ElementType  T;

        // view types from BS
        typedef typename BS::ConstView          ConstBSView;
        typedef typename BS::View               BSView;
        typedef typename BS::View               BSNoView;

        typedef typename BS::ConstVectorView    ConstBSVectorView;
        typedef typename BS::VectorView         BSVectorView;
        typedef typename BS::VectorNoView       BSVectorNoView;

        // view types for HbMatrix
        typedef DenseVector<ConstBSVectorView>  ConstVectorView;
        typedef DenseVector<BSVectorView>       VectorView;
        typedef DenseVector<BSVectorNoView>     VectorNoView;

        typedef GbMatrix<ConstBSView>           ConstGeneralView;
        typedef GbMatrix<BSView>                GeneralView;
        typedef GbMatrix<BSNoView>              GeneralNoView;

        typedef TbMatrix<ConstBSView>           ConstTriangularView;
        typedef TbMatrix<BSView>                TriangularView;
        typedef TbMatrix<BSNoView>              TriangularNoView;

        HbMatrix();

        HbMatrix(int dim, StorageUpLo upLo, int numOffDiags, int firstIndex=1);

        HbMatrix(const BS &bs, StorageUpLo upLo);

        HbMatrix(const HbMatrix<BS> &rhs);

        template <typename RHS>
            HbMatrix(const HbMatrix<RHS> &rhs);

        HbMatrix(const TbMatrix<BS> &rhs);

        template <typename RHS>
            HbMatrix(const TbMatrix<RHS> &rhs);

        // -- operators --------------------------------------------------------

        HbMatrix<BS> &
        operator*=(T alpha);

        HbMatrix<BS> &
        operator/=(T alpha);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        // -- views ------------------------------------------------------------
        ConstVectorView
        diag(int d) const;

        VectorView
        diag(int d);

        // general views
        ConstGeneralView
        general() const;

        GeneralView
        general();

        // triangular views
        ConstTriangularView
        triangular() const;

        TriangularView
        triangular();

        // -- methods ----------------------------------------------------------

        // for BLAS/LAPACK
        StorageUpLo
        upLo() const;

        int
        dim() const;

        int
        numOffDiags() const;

        int
        leadingDimension() const;

        const T *
        data() const;

        T *
        data();

        // for element access
        int
        firstIndex() const;

        int
        lastIndex() const;

        Range
        indices() const;

        Range
        diags() const;

        // -- implementation ---------------------------------------------------
        const BS &
        engine() const;

        BS &
        engine();

    private:
        BS _bs;
        StorageUpLo _upLo;
};

template <typename BS>
struct TypeInfo<HbMatrix<BS> >
{
    typedef HbMatrix<BS> Impl;
    typedef typename BS::ElementType ElementType;
};

// == HpMatrix =================================================================

template <typename PS>
class HpMatrix
    : public HermitianMatrix<HpMatrix<PS> >
{
    public:
        typedef typename HpMatrix<PS>::ElementType  T;

        // view types from PS
        typedef typename PS::ConstView          ConstPSView;
        typedef typename PS::View               PSView;
        typedef typename PS::View               PSNoView;

        // matrix-view types for TpMatrix
        typedef TpMatrix<ConstPSView>           ConstTriangularView;
        typedef TpMatrix<PSView>                TriangularView;
        typedef TpMatrix<PSNoView>              TriangularNoView;

        HpMatrix();

        HpMatrix(int dim, int firstIndex=1);

        HpMatrix(const PS &ps);

        // -- operators --------------------------------------------------------
        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        // -- views ------------------------------------------------------------

        // triangular views
        ConstTriangularView
        triangular(UnitDiag unitDiag=NonUnit) const;

        TriangularView
        triangular(UnitDiag unitDiag=NonUnit);

        // -- methods ----------------------------------------------------------

        // for BLAS/LAPACK
        StorageUpLo
        upLo() const;

        int
        dim() const;

        const T *
        data() const;

        T *
        data();

        // for element access
        int
        firstIndex() const;

        int
        lastIndex() const;

        Range
        indices() const;

        // -- implementation ---------------------------------------------------
        const PS &
        engine() const;

        PS &
        engine();

    private:
        PS _ps;
};

template <typename PS>
struct TypeInfo<HpMatrix<PS> >
{
    typedef HpMatrix<PS> Impl;
    typedef typename PS::ElementType ElementType;
};

} // namespace flens

#include <flens/hermitianmatrix.tcc>

#endif // FLENS_HERMITIANMATRIX_H
