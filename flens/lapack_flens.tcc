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

#include <cassert>
#include <complex>

#include <flens/lapack.h>
#include <flens/aux_complex.h>

namespace flens {

using std::complex;

//-- getrf ---------------------------------------------------------------------

template <typename FS>
int
trf(GeMatrix<FS> &A, DenseVector<Array<int> > &P)
{
    int info;
    info = getrf(A.numRows(), A.numCols(),
                 A.data(), A.leadingDimension(),
                 P.data());
    assert(info>=0);
    return info;
}

//-- gbtrf ---------------------------------------------------------------------

template <typename BS>
int
trf(GbMatrix<BS> &A, DenseVector<Array<int> > &P)
{
    int info;
    info = gbtrf(A.numRows(), A.numCols(),
                 A.numSubDiags(), A.numSuperDiags()-A.numSubDiags(),
                 A.data(), A.leadingDimension(),
                 P.data());
    assert(info>=0);
    return info;
}

//-- potrf ---------------------------------------------------------------------

template <typename FS>
int
trf(SyMatrix<FS> &A)
{
    return potrf(A.upLo(), A.dim(), A.data(),A.leadingDimension());
}

//-- getri ---------------------------------------------------------------------

template <typename MA>
int
tri(GeMatrix<MA> &A, DenseVector<Array<int> > &P)
{
    assert(A.numRows()==A.numCols());
    assert(P.length()==A.numCols());

    typedef typename MA::ElementType T;
    T lwork;
    int info;

    // query optimal work space size
    info = getri(A.numCols(), A.data(), A.leadingDimension(), P.data(),
                 &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(iabsround(lwork));
    return getri(A.numCols(), A.data(), A.leadingDimension(), P.data(),
                 work.data(), work.length());
}

//-- getrs ---------------------------------------------------------------------

template <typename MA, typename MB>
int
trs(Transpose trans, const GeMatrix<MA> &LU,
    const DenseVector<Array<int> > &P, GeMatrix<MB> &B)
{
    assert(LU.numRows()==LU.numCols());
    assert(P.length()==LU.numRows());

    return getrs(trans, LU.numCols(), B.numCols(),
                 LU.data(), LU.leadingDimension(),
                 P.data(), B.data(), B.leadingDimension());
}

template <typename MA, typename VB>
int
trs(Transpose trans, const GeMatrix<MA> &LU,
    const DenseVector<Array<int> > &P,
    DenseVector<VB> &B)
{
    assert(LU.numRows()==LU.numCols());
    assert(P.length()==LU.numRows());
    assert(B.stride()==1);

    return getrs(trans, LU.numCols(), 1,
                 LU.data(), LU.leadingDimension(),
                 P.data(), B.data(), B.length());
}

//-- gbtrs ---------------------------------------------------------------------

template <typename MA, typename MB>
int
trs(Transpose trans, GbMatrix<MA> &LU,
    const DenseVector<Array<int> > &P, GeMatrix<MB> &B)
{
    assert(LU.numRows()==LU.numCols());
    assert(P.length()==LU.numRows());

    return gbtrs(trans, LU.numRows(),
                 LU.numSubDiags(), LU.numSuperDiags()-LU.numSubDiags(),
                 B.numCols(), LU.data(), LU.leadingDimension(),
                 P.data(), B.data(), B.leadingDimension());
}

template <typename MA, typename VB>
int
trs(Transpose trans, const GbMatrix<MA> &LU,
    const DenseVector<Array<int> > &P, DenseVector<VB> &B)
{
    assert(LU.numRows()==LU.numCols());
    assert(P.length()==LU.numRows());
    assert(B.stride()==1);

    return gbtrs(trans, LU.numRows(),
                 LU.numSubDiags(), LU.numSuperDiags()-LU.numSubDiags(),
                 1, LU.data(), LU.leadingDimension(),
                 P.data(), B.data(), B.length());
}

//-- gesv ----------------------------------------------------------------------

template <typename MA, typename MB>
int
sv(GeMatrix<MA> &A, DenseVector<Array<int> > &P, GeMatrix<MB> &B)
{
    assert(A.numRows()==A.numCols());
    assert(B.numRows()==A.numRows());

    int info;
    info = gesv(A.numRows(), B.numCols(), A.data(), A.leadingDimension(),
                P.data(), B.data(), B.leadingDimension());
    assert(info>=0);
    return info;
}

template <typename MA, typename VB>
int
sv(GeMatrix<MA> &A, DenseVector<Array<int> > &P, DenseVector<VB> &B)
{
    assert(A.numRows()==A.numCols());
    assert(B.length()==A.numRows());
    assert(B.stride()==1);

    int info;
    info = gesv(A.numRows(), 1, A.data(), A.leadingDimension(),
                P.data(), B.data(), B.length());
    assert(info>=0);
return info;
}

//-- gbsv ----------------------------------------------------------------------

template <typename MA, typename MB>
int
sv(GbMatrix<MA> &A, DenseVector<Array<int> > &P, GeMatrix<MB> &B)
{
    assert(A.numRows()==A.numCols());

    return gbsv(A.numRows(),
                A.numSubDiags(), A.numSuperDiags()-A.numSubDiags(),
                B.numCols(),
                A.data(), A.leadingDimension(),
                P.data(), B.data(), B.leadingDimension());
}

template <typename MA, typename VB>
int
sv(GbMatrix<MA> &A, DenseVector<Array<int> > &P, DenseVector<VB> &B)
{
    assert(A.numRows()==A.numCols());

    return gbsv(A.numRows(),
                A.numSubDiags(), A.numSuperDiags()-A.numSubDiags(),
                1,
                A.data(), A.leadingDimension(),
                P.data(), B.data(), B.length());
}

//-- trtrs ---------------------------------------------------------------------

template <typename MA, typename MB>
int
trs(Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B)
{
    return trtrs(A.upLo(), trans, A.unitDiag(), A.dim(), B.numCols(),
                 A.data(), A.leadingDimension(),
                 B.data(), B.leadingDimension());
}

template <typename MA, typename VB>
int
trs(Transpose trans, const TrMatrix<MA> &A, DenseVector<VB> &B)
{
    assert(B.stride()==1);

    return trtrs(A.upLo(), trans, A.unitDiag(), A.dim(), 1,
                 A.data(), A.leadingDimension(), B.data(), B.length());
}

//-- geqrf ---------------------------------------------------------------------

template <typename MA, typename VT>
int
qrf(GeMatrix<MA> &A, DenseVector<VT> &tau)
{
    typedef typename MA::ElementType T;
    T lwork;
    int info;
    tau.resize(std::min(A.numRows(),A.numCols()));

    // query optimal work space size
    info = geqrf(A.numRows(), A.numCols(),
                 A.data(), A.leadingDimension(), tau.data(),
                 &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return geqrf(A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                 tau.data(), work.data(), work.length());
}

//-- orgqr ---------------------------------------------------------------------

template <typename MA, typename VT>
int
orgqr(GeMatrix<MA> &A, const DenseVector<VT> &tau)
{
    typedef typename MA::ElementType T;
    T lwork;
    int info;

    // query optimal work space size
    info = orgqr(A.numRows(), A.numCols(), tau.length(),
                 A.data(), A.leadingDimension(), tau.data(),
                 &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return orgqr(A.numRows(), A.numCols(), tau.length(),
                 A.data(), A.leadingDimension(), tau.data(),
                 work.data(), work.length());
}

//-- ormqr ---------------------------------------------------------------------

template <typename MA, typename VT, typename MC>
int
ormqr(BlasSide side, Transpose trans,
      const GeMatrix<MA> &A, const DenseVector<VT> &tau,
      GeMatrix<MC> &C)
{
    typedef typename MA::ElementType T;
    T lwork;
    int info;

    // query optimal work space size
    info = ormqr(side, trans, C.numRows(), C.numCols(), tau.length(),
                 A.data(), A.leadingDimension(), tau.data(),
                 C.data(), C.leadingDimension(),
                 &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return ormqr(side, trans, C.numRows(), C.numCols(), tau.length(),
                 A.data(), A.leadingDimension(), tau.data(),
                 C.data(), C.leadingDimension(),
                 work.data(), work.length());
}

//-- gels ----------------------------------------------------------------------

template <typename MA, typename MB>
int
ls(Transpose trans, GeMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename MA::ElementType T;
    T lwork;
    int info;

    // query work space size
    info = gels(trans, A.numRows(), A.numCols(), B.numCols(),
                A.data(), A.leadingDimension(),
                B.data(), B.leadingDimension(),
                &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return gels(trans, A.numRows(), A.numCols(), B.numCols(),
                A.data(), A.leadingDimension(),
                B.data(), B.leadingDimension(),
                work.data(), work.length());
}

//-- gelss ----------------------------------------------------------------------

template <typename MA, typename MB>
int
lss(GeMatrix<MA> &A, GeMatrix<MB> &B)
{
    typedef typename MA::ElementType T;
    T lwork;
    T rcond = -1.;
    int rank = -1;
    int info;

    // singular values output
    int ssize = A.numRows()<A.numCols() ? A.numRows() : A.numCols();
    DenseVector<Array<T> > s(ssize);

    // query work space size
    info = gelss(A.numRows(), A.numCols(), B.numCols(),
                 A.data(), A.leadingDimension(),
                 B.data(), B.leadingDimension(),
                 s.data(), rcond, rank, &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return gelss(A.numRows(), A.numCols(), B.numCols(),
                 A.data(), A.leadingDimension(),
                 B.data(), B.leadingDimension(),
                 s.data(), rcond, rank, work.data(), work.length());
}

//-- gees ----------------------------------------------------------------------

// real version

template <typename MA, typename WR, typename WI, typename MVS>
int
es(bool CalcEv, bool sort,
   typename ESSelectInfo<typename MA::ElementType>::pfunc select,
   GeMatrix<MA> &A, int &sdim,
   DenseVector<WR> &wr, DenseVector<WI> &wi,
   GeMatrix<MVS> &VS)
{
    typedef typename MA::ElementType T;
    int N = A.numCols();

    assert(A.numRows()==N);
    assert(wr.length()==N);
    assert(wi.length()==N);
    if (CalcEv) {
        assert (VS.numRows()==N);
        assert (VS.numCols()==N);
    }
    DenseVector< Array<int> > bwork;
    if (sort) {
        bwork.resize(N);
    }

    int * bwork_data = (sort ? bwork.data() : 0);
    T *    vs_data    = (CalcEv ? VS.data() : 0);
    int    ldvs       = (CalcEv ? VS.leadingDimension() : 1);

    // query optimal work space size
    T lwork;
    int info = gees(CalcEv, sort, select,
                    N, A.data(), A.leadingDimension(), sdim,
                    wr.data(), wi.data(), vs_data, ldvs,
                    &lwork, -1, bwork_data );
    assert (info==0);

    DenseVector< Array<T> > work(int(lwork.real()));
    return gees(CalcEv, sort, select,
                N, A.data(), A.leadingDimension(), sdim,
                wr.data(), wi.data(), vs_data, ldvs,
                work.data(), work.length(), bwork_data );
}

// complex version

template <typename MA, typename VW, typename MVS>
int
es(bool CalcEv, bool sort,
   typename ESSelectInfo<typename MA::ElementType>::pfunc select,
   GeMatrix<MA> &A, int &sdim,
   DenseVector<VW> &w, GeMatrix<MVS> &VS)
{
    typedef typename MA::ElementType T;
    typedef typename MA::ElementType::value_type RT;
    int N = A.numCols();

    assert(A.numRows()==N);
    assert(w.length()==N);
    if (CalcEv) {
        assert (VS.numRows()==N);
        assert (VS.numCols()==N);
    }
    DenseVector< Array<int> > bwork;
    if (sort) {
        bwork.resize(N);
    }

    int * bwork_data = (sort ? bwork.data() : 0);
    T *    vs_data    = (CalcEv ? VS.data() : 0);
    int    ldvs       = (CalcEv ? VS.leadingDimension() : 1);

    DenseVector< Array<RT> > rwork(A.numCols());
    // query optimal work space size
    T lwork;
    int info = gees(CalcEv, sort, select,
                    N, A.data(), A.leadingDimension(), sdim, w.data(),
                    vs_data, ldvs,
                    &lwork, -1, rwork.data(), bwork_data );
    assert (info==0);

    DenseVector< Array<T> > work(int(lwork.real()));
    return gees(CalcEv, sort, select,
                N, A.data(), A.leadingDimension(), sdim, w.data(),
                vs_data, ldvs,
                work.data(), work.length(), rwork.data(), bwork_data );
}

//-- geev ----------------------------------------------------------------------

template <typename MA, typename WR, typename WI, typename VL, typename VR>
int
ev(bool leftEV, bool rightEV,
   GeMatrix<MA> &A, DenseVector<WR> &wr, DenseVector<WI> &wi,
   GeMatrix<VL> &vl, GeMatrix<VR> &vr)
{
    assert(A.numRows()==A.numCols());
    assert(wr.length()==A.numRows());
    assert(wi.length()==A.numRows());
    assert(vl.numRows()==vl.numCols());
    assert(vr.numRows()==vr.numCols());
    assert(!leftEV || (vl.numRows()==A.numRows()));
    assert(!rightEV || (vr.numRows()==A.numRows()));

    typedef typename MA::ElementType T;
    T lwork;
    int info;

    int ldvl = leftEV  ? vl.engine().leadingDimension() : 1;
    int ldvr = rightEV ? vr.engine().leadingDimension() : 1;
    T *vldata = leftEV  ? vl.data() : 0;
    T *vrdata = rightEV ? vr.data() : 0;

    // query optimal work space size
    info = geev(leftEV, rightEV, A.numRows(), A.data(), A.leadingDimension(),
                wr.data(), wi.data(),
                vldata, ldvl,
                vrdata, ldvr,
                &lwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork));
    return geev(leftEV, rightEV, A.numRows(), A.data(), A.leadingDimension(),
                wr.data(), wi.data(),
                vldata, ldvl,
                vrdata, ldvr,
                work.data(), work.length());
}

template <typename MA, typename W, typename VL, typename VR>
int
ev(bool leftEV, bool rightEV,
   GeMatrix<MA> &A, DenseVector<W> &w, GeMatrix<VL> &vl, GeMatrix<VR> &vr)
{
    assert(A.numRows()==A.numCols());
    assert(w.length()==A.numRows());
    assert(vl.numRows()==vl.numCols());
    assert(vr.numRows()==vr.numCols());
    assert(!leftEV || (vl.numRows()==A.numRows()));
    assert(!rightEV || (vr.numRows()==A.numRows()));

    int ldvl = leftEV  ? vl.engine().leadingDimension() : 1;
    int ldvr = rightEV ? vr.engine().leadingDimension() : 1;
    typename VL::ElementType *vldata = leftEV  ? vl.data() : 0;
    typename VR::ElementType *vrdata = rightEV ? vr.data() : 0;
    typedef typename MA::ElementType T;
    T lwork;
    DenseVector<Array<typename T::value_type> > rwork(2*A.numRows());
    int info;

    // query optimal work space size
    info = geev(leftEV, rightEV, A.numRows(), A.data(), A.leadingDimension(),
                w.data(),
                vldata, ldvl,
                vrdata, ldvr,
                &lwork, -1, rwork.data());
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(lwork.real()));
    return geev(leftEV, rightEV, A.numRows(), A.data(), A.leadingDimension(),
                w.data(),
                vldata, ldvl,
                vrdata, ldvr,
                work.data(), work.length(), rwork.data());
}

//-- syev ----------------------------------------------------------------------

template <typename MA, typename VW>
int
ev(bool computeEigenvectors, SyMatrix<MA> &A, DenseVector<VW> &w)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }

    typedef typename MA::ElementType T;
    T optwork;
    int info;

    // query optimal work space size
    info = syev(computeEigenvectors, A.upLo(), A.dim(),
                A.data(), A.leadingDimension(), w.data(),
                &optwork, -1);
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(static_cast<int>(optwork));
    return syev(computeEigenvectors, A.upLo(), A.dim(),
                A.data(), A.leadingDimension(), w.data(),
                work.data(), work.length());
}

//-- sbev ----------------------------------------------------------------------

template <typename MA, typename VW, typename MZ>
int
ev(bool computeEigenvectors, SbMatrix<MA> &A, DenseVector<VW> &w,
   GeMatrix<MZ> &Z)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }
    if (computeEigenvectors) {
        if ((Z.numRows()!=A.dim()) || (Z.numCols()!=A.dim())) {
            Z.resize(A.dim(),A.dim());
        }
    }

    typedef typename MA::ElementType T;
    DenseVector<Array<T> > work(3*A.dim()-2);
    return sbev(computeEigenvectors,
                A.upLo(), A.dim(), A.numOffDiags(),
                A.data(), A.leadingDimension(),
                w.data(), Z.data(), Z.leadingDimension(), work.data());
}

//-- spev ----------------------------------------------------------------------

template <typename MA, typename VW, typename MZ>
int
ev(bool computeEigenvectors, SpMatrix<MA> &A, DenseVector<VW> &w,
   GeMatrix<MZ> &Z)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }
    if (computeEigenvectors) {
        if ((Z.numRows()!=A.dim()) || (Z.numCols()!=A.dim())) {
            Z.resize(A.dim(),A.dim());
        }
    }
    typedef typename MA::ElementType T;
    DenseVector<Array<T> > work(3*A.dim());
    return spev(computeEigenvectors,
                A.upLo(), A.dim(), A.data(), w.data(),
                Z.data(), Z.leadingDimension(), work.data());
}

//-- heev ----------------------------------------------------------------------

template <typename MA, typename VW>
int
ev(bool computeEigenvectors, HeMatrix<MA> &A, DenseVector<VW> &w)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }

    typedef typename MA::ElementType T;
    T optwork;
    DenseVector<Array<typename T::value_type> > rwork(3*A.dim()-2);
    int info;

    // query optimal work space size
    info = heev(computeEigenvectors, A.upLo(), A.dim(),
                A.data(), A.leadingDimension(),
                w.data(),
                &optwork, -1,
                rwork.data());
    assert(info==0);

    // allocate work space
    DenseVector<Array<T> > work(int(optwork.real()));
    return heev(computeEigenvectors, A.upLo(), A.dim(),
                A.data(), A.leadingDimension(),
                w.data(),
                work.data(), work.length(),
                rwork.data());
}

//-- hbev ----------------------------------------------------------------------

template <typename MA, typename VW, typename MZ>
int
ev(bool computeEigenvectors, HbMatrix<MA> &A, DenseVector<VW> &w,
   GeMatrix<MZ> &Z)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }
    if (computeEigenvectors) {
        if ((Z.numRows()!=A.dim()) || (Z.numCols()!=A.dim())) {
            Z.resize(A.dim(),A.dim());
        }
    }

    typedef typename MA::ElementType T;
    DenseVector<Array<T> > work(A.dim());
    DenseVector<Array<typename T::value_type> > rwork(3*A.dim()-2);
    return hbev(computeEigenvectors,
                A.upLo(), A.dim(), A.numOffDiags(),
                A.data(), A.leadingDimension(),
                w.data(), Z.data(), Z.leadingDimension(),
                work.data(), rwork.data());
}

//-- hpev ----------------------------------------------------------------------

template <typename MA, typename VW, typename MZ>
int
ev(bool computeEigenvectors, HpMatrix<MA> &A, DenseVector<VW> &w,
   GeMatrix<MZ> &Z)
{
    if (w.length()!=A.dim()) {
        w.resize(A.dim());
    }
    if (computeEigenvectors) {
        if ((Z.numRows()!=A.dim()) || (Z.numCols()!=A.dim())) {
            Z.resize(A.dim(),A.dim());
        }
    }

    typedef typename MA::ElementType T;
    DenseVector<Array<T> > work(2*A.dim()-1);
    DenseVector<Array<typename T::value_type> > rwork(3*A.dim()-2);
    return hpev(computeEigenvectors, A.upLo(), A.dim(), A.data(),
                w.data(), Z.data(), Z.leadingDimension(),
                work.data(), rwork.data());
}

//-- gesvd ---------------------------------------------------------------------
template <typename MA, typename VS, typename MU, typename MV>
int
svd(GeMatrix<MA> &A,
    DenseVector<VS> &s,
    GeMatrix<MU> &U,
    GeMatrix<MV> &V)
{
    return svd(All,All,A,s,U,V);
}

// Template hack to support implementations depending on element type
// (c++ does not allow partial specialization of func templates)

// non-complex implementation

template <typename T, typename MA, typename VS, typename VU, typename VVT>
struct SvdImpl {
    static int
    func(SVectorsJob jobu, SVectorsJob jobvt,
         GeMatrix<MA> &A,
         DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
    {
        int m = A.numRows();
        int n = A.numCols();
        int mindim = (m<n ? m : n);
        if (S.length()!=mindim) {
            S.resize(mindim);
        }

        switch (jobu) {
            case All:
                if ((U.numRows()!=m) || (U.numCols()!=m)) {
                    U.resize(m,m);
                }
                break;
            case SmallDim:
                if ((U.numRows()!=m) || (U.numCols()!=mindim)) {
                    U.resize(m,mindim);
                }
                break;
            case Overwrite:
                assert(0); // case 'Overwrite' not allowed!
            case None:
                assert(0); // case 'None' not allowed!
        }

        switch (jobvt) {
            case All:
                if ((VT.numRows()!=n) || (VT.numCols()!=n)) {
                    VT.resize(n,n);
                }
                break;
            case SmallDim:
                if ((VT.numRows()!=mindim) || (VT.numCols()!=n)) {
                    VT.resize(mindim,n);
                }
                break;
            case Overwrite:
                assert(0); // case 'Overwrite' not allowed!
            case None:
                assert(0); // case 'None' not allowed!
        }

        // query work space size
        T lwork;
        int
        info = gesvd(jobu, jobvt,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     &lwork, -1 );
        assert(info==0);

        // allocate work space
        DenseVector<Array<T> > work(static_cast<int>(lwork));
        return gesvd(jobu, jobvt,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     work.data(), work.length() );
    }
};

// complex specialization

template <typename T, typename MA, typename VS, typename VU, typename VVT>
struct SvdImpl<complex<T>, MA,VS,VU,VVT> {
    static int
    func(SVectorsJob jobu, SVectorsJob jobvt,
         GeMatrix<MA> &A,
         DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
    {
        int m = A.numRows();
        int n = A.numCols();
        int mindim = (m<n ? m : n);
        if (S.length()!=mindim) {
            S.resize(mindim);
        }
        assert(jobu!=Overwrite || jobvt!=Overwrite);

        switch (jobu) {
            case All:
                if ((U.numRows()!=m) || (U.numCols()!=m)) {
                    U.resize(m,m);
                }
                break;
            case SmallDim:
                if ((U.numRows()!=m) || (U.numCols()!=mindim)) {
                    U.resize(m,mindim);
                }
                break;
            case Overwrite:
                assert(0); // case 'Overwrite' not allowed!
            case None:
                assert(0); // case 'None' not allowed!
        }

        switch (jobvt) {
            case All:
                if ((VT.numRows()!=n) || (VT.numCols()!=n)) {
                    VT.resize(n,n);
                }
                break;
            case SmallDim:
                if ((VT.numRows()!=mindim) || (VT.numCols()!=n)) {
                    VT.resize(mindim,n);
                }
                break;
            case Overwrite:
                assert(0); // case 'Overwrite' not allowed!
            case None:
                assert(0); // case 'None' not allowed!
        }

        DenseVector< Array<T> > rwork(5*mindim);

        // query work space size
        complex<T> lwork;
        int
        info = gesvd(jobu, jobvt,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     &lwork, -1, rwork.data() );
        assert(info==0);

        // allocate work space
        DenseVector< Array<complex<T> > > work(static_cast<int>(lwork.real()));
        return gesvd(jobu, jobvt,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     work.data(), work.length(), rwork.data() );
    }
};

template <typename MA, typename VS, typename VU, typename VVT>
int
svd(SVectorsJob jobu, SVectorsJob jobvt,
    GeMatrix<MA> &A, DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
{
    typedef SvdImpl<typename MA::ElementType, MA,VS,VU,VVT> impl;
    return impl::func(jobu, jobvt, A, S, U, VT);
}

//-- gesdd ---------------------------------------------------------------------

// Template hack to support implementations depending on element type
// (c++ does not allow partial specialization of func templates)

// non-complex implementation

template <typename T, typename MA, typename VS, typename VU, typename VVT>
struct SddImpl {
    static int
    func(SVectorsJob jobz,
         GeMatrix<MA> &A,
         DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
    {
        int mindim = (A.numRows()<A.numCols() ? A.numRows() : A.numCols());
        assert(S.length()==mindim);

        switch (jobz) {
            case All:
                assert(U.numRows()==A.numRows());
                assert(U.numCols()==A.numRows());
                assert(VT.numRows()==A.numCols());
                assert(VT.numCols()==A.numCols());
                break;
            case SmallDim:
                assert(U.numRows()==A.numRows());
                assert(U.numCols()==mindim);
                assert(VT.numRows()==mindim);
                assert(VT.numCols()==A.numCols);
                break;
            case Overwrite:
                if (mindim<A.numCols()) {
                    assert(U.numRows()==A.numRows());
                    assert(U.numCols()==A.numRows());
                } else {
                    assert(VT.numRows()==A.numcols());
                    assert(VT.numCols()==A.numcols());
                }
                break;
            case None:
                assert(0); // case 'None' not allowed!
        }

        DenseVector< Array<int> > iwork(8*mindim);

        // query work space size
        T lwork;
        int
        info = gesdd(jobz,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     &lwork, -1 , iwork.data() );
        assert(info==0);

        // allocate work space
        DenseVector<Array<T> > work((int) lwork);
        return gesvd(jobz,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     work.data(), work.length(), iwork);
    }
};

// complex specialization

template <typename T, typename MA, typename VS, typename VU, typename VVT>
struct SddImpl<complex<T>, MA,VS,VU,VVT> {
    static int
    func(SVectorsJob jobz,
         GeMatrix<MA> &A,
         DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
    {
        int mindim = (A.numRows()<A.numCols() ? A.numRows() : A.numCols());
        assert(S.length()==mindim);

        int lrwork = 5*(mindim*mindim + mindim);
        switch (jobz) {
            case All:
                assert(U.numRows()==A.numRows());
                assert(U.numCols()==A.numRows());
                assert(VT.numRows()==A.numCols());
                assert(VT.numCols()==A.numCols());
                break;
            case SmallDim:
                assert(U.numRows()==A.numRows());
                assert(U.numCols()==mindim);
                assert(VT.numRows()==mindim);
                assert(VT.numCols()==A.numCols());
                break;
            case Overwrite:
                if (mindim<A.numCols()) {
                    assert(U.numRows()==A.numRows());
                    assert(U.numCols()==A.numRows());
                } else {
                    assert(VT.numRows()==A.numCols());
                    assert(VT.numCols()==A.numCols());
                }
                break;
            case None:
                lrwork = 7*mindim;
        }

        DenseVector< Array<T> > rwork(lrwork);
        DenseVector< Array<int> > iwork(8*mindim);

        // query work space size
        complex<T> lwork;
        int
        info = gesdd(jobz,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     &lwork, -1, rwork.data(), iwork.data() );
        assert(info==0);

        // allocate work space
        DenseVector< Array<complex<T> > > work((int) lwork.real());
        return gesdd(jobz,
                     A.numRows(), A.numCols(), A.data(), A.leadingDimension(),
                     S.data(),
                     U.data(), U.leadingDimension(),
                     VT.data(), VT.leadingDimension(),
                     work.data(), work.length(), rwork.data(), iwork.data() );
    }
};

template <typename MA, typename VS, typename VU, typename VVT>
int
sdd(SVectorsJob jobz,
    GeMatrix<MA> &A, DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT)
{
    typedef SddImpl<typename MA::ElementType, MA,VS,VU,VVT> impl;
    return impl::func(jobz, A, S, U, VT);
}

} // namespace flens
