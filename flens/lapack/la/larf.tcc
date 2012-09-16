/*
 *   Copyright (c) 2011, Michael Lehn
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

/* Based on
 *     SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LARF_TCC
#define FLENS_LAPACK_LA_LARF_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larf_impl(Side side, const DenseVector<VV> &v, const TAU &tau,
          GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    using lapack::ilalc;
    using lapack::ilalr;

    typedef typename GeMatrix<MC>::ElementType  ElementType;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const ElementType  Zero(0), One(1);

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();

    IndexType lastV = 0;
    IndexType lastC = 0;

    if (tau!=Zero) {
//
//      Set up variables for scanning V.  LASTV begins pointing to the end of V.
//      Look for the last non-zero row in V.
//
        if (side==Left) {
            lastV = m;
        } else {
            lastV = n;
        }
        IndexType i = (v.inc()>0) ? 1 + (lastV-1)*v.inc() : 1;
//
//      Look for the last non-zero row in V.
//
        while (lastV>0 && v(i)==Zero) {
            --lastV;
            i -= v.inc();
        }
        if (side==Left) {
//
//          Scan for the last non-zero column in C(1:lastv,:).
//
            lastC = ilalc(C(_(1,lastV),_));
        } else {
//
//          Scan for the last non-zero row in C(:,1:lastv).
//
            lastC = ilalr(C(_,_(1,lastV)));
        }
    }
//
//  Note that lastc.eq.0 renders the BLAS operations null; no special
//  case is needed at this level.
//
    const auto _v = v(_(1,lastV));

    if (side==Left) {
//
//      Form  H * C
//
        if (lastV>0) {
            auto _work = work(_(1,lastC));
            auto _C = C(_(1,lastV),_(1,lastC));
//
//          work(1:lastc,1) := C(1:lastv,1:lastc)' * v(1:lastv,1)
//
            blas::mv(Trans, One, _C, _v, Zero, _work);
//
//          C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * work(1:lastc,1)'
//
            blas::r(-tau, _v, _work, _C);
        }
    } else {
//
//      Form  C * H
//
        if (lastV>0) {
            auto _work = work(_(1,lastC));
            auto _C = C(_(1,lastC),_(1,lastV));
//
//          work(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
//
            blas::mv(NoTrans, One, _C, _v, Zero, _work);
//
//          C(1:lastc,1:lastv) := C(...) - work(1:lastc,1) * v(1:lastv,1)'
//
            blas::r(-tau, _work, _v, _C);
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larf_impl(Side side, const DenseVector<VV> &v, const TAU &tau,
          GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    typedef typename GeMatrix<MC>::IndexType  IndexType;

    cxxlapack::larf<IndexType>(getF77Char(side),
                               C.numRows(),
                               C.numCols(),
                               v.data(),
                               v.inc()*v.stride(),
                               tau,
                               C.data(),
                               C.leadingDimension(),
                               work.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larf(Side side, const DenseVector<VV> &v, const TAU &tau,
     GeMatrix<MC> &C, DenseVector<VWORK> &work)
{
    LAPACK_DEBUG_OUT("larf");

//
//  Test the input parameters
//
    ASSERT((v.inc()>0 && v.firstIndex()==1)
        || (v.inc()<0 && v.lastIndex()==1));
    ASSERT((side==Left && v.length()==C.numRows())
        || (side==Right && v.length()==C.numCols()));
    ASSERT(C.firstRow()==1);
    ASSERT(C.firstCol()==1);
    ASSERT((side==Left && work.length()>=C.numCols())
        || (side==Right && work.length()>=C.numRows()));

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MC>::NoView       C_org      = C;
    typename DenseVector<VV>::NoView    work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::larf_impl(side, v, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MC>::NoView       C_generic      = C;
    typename DenseVector<VV>::NoView    work_generic   = work;

//
//  restore output arguments
//
    C       = C_org;
    work    = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::larf_impl(side, v, tau, C, work);

    bool failed = false;

    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename VV, typename TAU, typename MC, typename VWORK>
void
larf(Side side, const VV &v, const TAU &tau, MC &&C, VWORK &&work)
{
    CHECKPOINT_ENTER;
    larf(side, v, tau, C, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LARF_TCC
