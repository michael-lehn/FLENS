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

namespace flens {

// == SyMatrix =================================================================

template <typename FS>
HeMatrix<FS>::HeMatrix()
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(int dim, StorageUpLo upLo, int firstIndex)
    : _fs(dim, dim, firstIndex, firstIndex), _upLo(upLo)
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const FS &fs, StorageUpLo upLo)
    : _fs(fs), _upLo(upLo)
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const HeMatrix<FS> &rhs)
    : HermitianMatrix<HeMatrix<FS> >(),
      _fs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const HeMatrix<RHS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename FS>
HeMatrix<FS>::HeMatrix(const TrMatrix<FS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

template <typename FS>
template <typename RHS>
HeMatrix<FS>::HeMatrix(const TrMatrix<RHS> &rhs)
    : _fs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

// -- operators ----------------------------------------------------------------

template <typename FS>
HeMatrix<FS> &
HeMatrix<FS>::operator*=(T alpha)
{
    scal(alpha, *this);
    return *this;
}

template <typename FS>
HeMatrix<FS> &
HeMatrix<FS>::operator/=(T alpha)
{
    scal(T(1)/alpha, *this);
    return *this;
}

template <typename FS>
const typename HeMatrix<FS>::T &
HeMatrix<FS>::operator()(int row, int col) const
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
    } else {
        assert(col<=row);
    }
#endif
    return _fs(row, col);
}

template <typename FS>
typename HeMatrix<FS>::T &
HeMatrix<FS>::operator()(int row, int col)
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
    } else {
        assert(col<=row);
    }
#endif
    return _fs(row, col);
}

// -- views --------------------------------------------------------------------

// general views
template <typename FS>
typename HeMatrix<FS>::ConstGeneralView
HeMatrix<FS>::general() const
{
    return ConstGeneralView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                            _fs.lastCol(), _fs.lastRow(),
                            firstRow(), firstCol()));
}

template <typename FS>
typename HeMatrix<FS>::GeneralView
HeMatrix<FS>::general()
{
    return GeneralView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                       _fs.lastCol(), _fs.lastRow(),
                       firstRow(), firstCol()));
}

// triangular views
template <typename FS>
typename HeMatrix<FS>::ConstTriangularView
HeMatrix<FS>::triangular() const
{
    return ConstTriangularView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                               _fs.lastCol(), _fs.lastRow(),
                               firstRow(), firstCol()), upLo(), NonUnit);
}

template <typename FS>
typename HeMatrix<FS>::TriangularView
HeMatrix<FS>::triangular()
{
    return TriangularView(_fs.view(_fs.firstRow(), _fs.firstCol(),
                          _fs.lastCol(), _fs.lastRow(),
                          firstRow(), firstCol()), upLo(), NonUnit);
}

// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename FS>
StorageUpLo
HeMatrix<FS>::upLo() const
{
    return _upLo;
}

template <typename FS>
int
HeMatrix<FS>::dim() const
{
    assert(_fs.numRows()==_fs.numCols());

    return _fs.numRows();
}

template <typename FS>
int
HeMatrix<FS>::leadingDimension() const
{
    return _fs.leadingDimension();
}

template <typename FS>
const typename HeMatrix<FS>::T *
HeMatrix<FS>::data() const
{
    return _fs.data();
}

template <typename FS>
typename HeMatrix<FS>::T *
HeMatrix<FS>::data()
{
    return _fs.data();
}

// for element access
template <typename FS>
int
HeMatrix<FS>::firstRow() const
{
    return _fs.firstRow();
}

template <typename FS>
int
HeMatrix<FS>::lastRow() const
{
    return _fs.lastRow();
}

template <typename FS>
int
HeMatrix<FS>::firstCol() const
{
    return _fs.firstCol();
}

template <typename FS>
int
HeMatrix<FS>::lastCol() const
{
    return _fs.lastCol();
}

template <typename FS>
Range
HeMatrix<FS>::rows()  const
{
    return _(firstRow(), lastRow());
}

template <typename FS>
Range
HeMatrix<FS>::cols() const
{
    return _(firstCol(), lastCol());
}

// -- implementation -----------------------------------------------------------

template <typename FS>
const FS &
HeMatrix<FS>::engine() const
{
    return _fs;
}

template <typename FS>
FS &
HeMatrix<FS>::engine()
{
    return _fs;
}

// == HbMatrix =================================================================

template <typename BS>
HbMatrix<BS>::HbMatrix()
{
}

template <typename BS>
HbMatrix<BS>::HbMatrix(int dim, StorageUpLo upLo, int numOffDiags,
                       int firstIndex)
    : _bs(dim, dim,
          (upLo==Lower) ? numOffDiags : 0,
          (upLo==Upper) ? numOffDiags : 0,
          firstIndex),
      _upLo(upLo)
{
}

template <typename BS>
HbMatrix<BS>::HbMatrix(const BS &bs, StorageUpLo upLo)
    : _bs(bs), _upLo(upLo)
{
}

template <typename BS>
HbMatrix<BS>::HbMatrix(const HbMatrix<BS> &rhs)
    : HermitianMatrix<HbMatrix<BS> >(),
      _bs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename BS>
template <typename RHS>
HbMatrix<BS>::HbMatrix(const HbMatrix<RHS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
}

template <typename BS>
HbMatrix<BS>::HbMatrix(const TbMatrix<BS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

template <typename BS>
template <typename RHS>
HbMatrix<BS>::HbMatrix(const TbMatrix<RHS> &rhs)
    : _bs(rhs.engine()), _upLo(rhs.upLo())
{
    assert(rhs.unitDiag()==NonUnit);
}

// -- operators ----------------------------------------------------------------

template <typename BS>
HbMatrix<BS> &
HbMatrix<BS>::operator*=(T alpha)
{
    scal(alpha, *this);
    return *this;
}

template <typename BS>
HbMatrix<BS> &
HbMatrix<BS>::operator/=(T alpha)
{
    scal(T(1)/alpha, *this);
    return *this;
}

template <typename BS>
const typename HbMatrix<BS>::T &
HbMatrix<BS>::operator()(int row, int col) const
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
        assert(col-row<=numOffDiags());
    } else {
        assert(col<=row);
        assert(row-col<=numOffDiags());
    }
#endif
    return _bs(row, col);
}

template <typename BS>
typename HbMatrix<BS>::T &
HbMatrix<BS>::operator()(int row, int col)
{
#ifndef NDEBUG
    if (_upLo==Upper) {
        assert(col>=row);
        assert(col-row<=numOffDiags());
    } else {
        assert(col<=row);
        assert(row-col<=numOffDiags());
    }
#endif
    return _bs(row, col);
}

// -- views --------------------------------------------------------------------

template <typename BS>
typename HbMatrix<BS>::ConstVectorView
HbMatrix<BS>::diag(int d) const
{
    return _bs.viewDiag(d);
}

template <typename BS>
typename HbMatrix<BS>::VectorView
HbMatrix<BS>::diag(int d)
{
    return _bs.viewDiag(d);
}

// general views
template <typename BS>
typename HbMatrix<BS>::ConstGeneralView
HbMatrix<BS>::general() const
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return ConstGeneralView(view);
}

template <typename BS>
typename HbMatrix<BS>::GeneralView
HbMatrix<BS>::general()
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return GeneralView(view);
}

// triangular views
template <typename BS>
typename HbMatrix<BS>::ConstTriangularView
HbMatrix<BS>::triangular() const
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return ConstTriangularView(view, upLo());
}

template <typename BS>
typename HbMatrix<BS>::TriangularView
HbMatrix<BS>::triangular()
{
    int firstDiag = -engine().numSubDiags();
    int lastDiag = engine().numSuperDiags();
    BSView view = engine().viewDiags(firstDiag, lastDiag, firstIndex());
    return TriangularView(view, upLo());
}

// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename BS>
StorageUpLo
HbMatrix<BS>::upLo() const
{
    return _upLo;
}

template <typename BS>
int
HbMatrix<BS>::dim() const
{
    assert (_bs.numRows()==_bs.numCols());

    return _bs.numRows();
}

template <typename BS>
int
HbMatrix<BS>::numOffDiags() const
{
    return (_upLo==Upper) ? _bs.numSuperDiags()
                          : _bs.numSubDiags();
}

template <typename BS>
int
HbMatrix<BS>::leadingDimension() const
{
    return _bs.leadingDimension();
}

template <typename BS>
const typename HbMatrix<BS>::T *
HbMatrix<BS>::data() const
{
    return _bs.data();
}

template <typename BS>
typename HbMatrix<BS>::T *
HbMatrix<BS>::data()
{
    return _bs.data();
}

// for element access
template <typename BS>
int
HbMatrix<BS>::firstIndex() const
{
    return _bs.firstRow();
}

template <typename BS>
int
HbMatrix<BS>::lastIndex() const
{
    assert(_bs.lastRow()==_bs.lastCol());
    return _bs.lastRow();
}

template <typename BS>
Range
HbMatrix<BS>::indices()  const
{
    return _(firstIndex(), lastIndex());
}

template <typename BS>
Range
HbMatrix<BS>::diags() const
{
    return (_upLo==Upper) ? _(0, numOffDiags())
                          : _(-numOffDiags(),0);
}

// -- implementation -----------------------------------------------------------

template <typename BS>
const BS &
HbMatrix<BS>::engine() const
{
    return _bs;
}

template <typename BS>
BS &
HbMatrix<BS>::engine()
{
    return _bs;
}

// == HpMatrix =================================================================

template <typename PS>
HpMatrix<PS>::HpMatrix()
{
}

template <typename PS>
HpMatrix<PS>::HpMatrix(int dim, int firstIndex)
    : _ps(dim, firstIndex)
{
}

template <typename PS>
HpMatrix<PS>::HpMatrix(const PS &ps)
    : _ps(ps)
{
}

// -- operators ----------------------------------------------------------------

template <typename PS>
const typename HpMatrix<PS>::T &
HpMatrix<PS>::operator()(int row, int col) const
{
    return _ps(row, col);
}

template <typename PS>
typename HpMatrix<PS>::T &
HpMatrix<PS>::operator()(int row, int col)
{
    return _ps(row, col);
}

// -- views --------------------------------------------------------------------

// triangular views
template <typename PS>
typename HpMatrix<PS>::ConstTriangularView
HpMatrix<PS>::triangular(UnitDiag unitDiag) const
{
    return ConstTriangularView(_ps.view(), unitDiag);
}

template <typename PS>
typename HpMatrix<PS>::TriangularView
HpMatrix<PS>::triangular(UnitDiag unitDiag)
{
    return TriangularView(_ps.view(), unitDiag);
}

// -- methods ------------------------------------------------------------------

// for BLAS/LAPACK
template <typename PS>
StorageUpLo
HpMatrix<PS>::upLo() const
{
    return StorageInfo<PS>::upLo;
}

template <typename PS>
int
HpMatrix<PS>::dim() const
{
    return _ps.dim();
}

template <typename PS>
const typename HpMatrix<PS>::T *
HpMatrix<PS>::data() const
{
    return _ps.data();
}

template <typename PS>
typename HpMatrix<PS>::T *
HpMatrix<PS>::data()
{
    return _ps.data();
}

// for element access
template <typename PS>
int
HpMatrix<PS>::firstIndex() const
{
    return _ps.firstIndex();
}

template <typename PS>
int
HpMatrix<PS>::lastIndex() const
{
    return _ps.lastIndex();
}

template <typename PS>
Range
HpMatrix<PS>::indices() const
{
    return _(firstIndex(), lastIndex());
}

// -- implementation -----------------------------------------------------------

template <typename PS>
const PS &
HpMatrix<PS>::engine() const
{
    return _ps;
}

template <typename PS>
PS &
HpMatrix<PS>::engine()
{
    return _ps;
}

} // namespace flens
