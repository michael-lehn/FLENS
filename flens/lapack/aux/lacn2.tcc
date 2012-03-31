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
 *
      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 *
 */

#ifndef FLENS_LAPACK_AUX_LACN2_TCC
#define FLENS_LAPACK_AUX_LACN2_TCC 1

#include <cmath>
#include <limits>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
template <typename  VV, typename VX, typename VSGN, typename EST,
          typename KASE, typename VSAVE>
void
lacn2_generic(DenseVector<VV> &v, DenseVector<VX> &x, DenseVector<VSGN> &sgn,
              EST &est, KASE &kase, DenseVector<VSAVE> &iSave)
{
    using std::abs;

    typedef typename DenseVector<VV>::ElementType   T;
    typedef typename DenseVector<VV>::IndexType     IndexType;
    const T Zero(0), One(1), Two(2);

    const IndexType itMax = 5;
    const IndexType n = v.length();
//  ..
//  .. Executable Statements ..
//
    if (kase==0) {
        for (IndexType i=1; i<=n; ++i) {
            x(i) = One / T(n);
        }
        kase = 1;
        iSave(1) = 1;
        return;
    }

    IndexType jLast;
    T         altSign, estOld;
    bool      converged;

    switch (iSave(1)) {
//
//  ................ ENTRY   (ISAVE( 1 ) = 1)
//  FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
//
    case 1:
        if (n==1) {
            v(1) = x(1);
            est = abs(v(1));
    //      ... QUIT
            goto QUIT;
        }
        est = blas::asum(x);

        for (IndexType i=1; i<=n; ++i) {
            x(i)    = sign(One, x(i));
            sgn(i)  = nint(x(i));
        }
        kase = 2;
        iSave(1) = 2;
        return;
//
//  ................ ENTRY   (ISAVE( 1 ) = 2)
//  FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
//
    case 2:
        iSave(2) = blas::iamax(x);
        iSave(3) = 2;
//
//      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
//
    MAIN:
        for (IndexType i=1; i<=n; ++i) {
            x(i) = Zero;
        }
        x(iSave(2)) = One;
        kase = 1;
        iSave(1) = 3;
        return;
//
//  ................ ENTRY   (ISAVE( 1 ) = 3)
//  X HAS BEEN OVERWRITTEN BY A*X.
//
    case 3:
        v = x;
        estOld = est;
        est = blas::asum(v);

        converged = true;
        for (IndexType i=1; i<=n; ++i) {
            if (nint(sign(One,x(i)))!=sgn(i)) {
                converged = false;
                break;
            }
        }
        if (converged) {
//          REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
            goto ITERATION_COMPLETE;
        }

//      TEST FOR CYCLING.
        if (est<=estOld) {
            goto ITERATION_COMPLETE;
        }

        for (IndexType i=1; i<=n; ++i) {
            x(i)    = sign(One, x(i));
            sgn(i)  = nint(x(i));
        }
        kase = 2;
        iSave(1) = 4;
        return;
//
//  ................ ENTRY   (ISAVE( 1 ) = 4)
//  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
//
    case 4:
        jLast = iSave(2);
        iSave(2) = blas::iamax(x);
        if (x(jLast)!=abs(x(iSave(2))) && iSave(3)<itMax) {
            ++iSave(3);
            goto MAIN;
        }
//
//      ITERATION COMPLETE.  FINAL STAGE.
//
    ITERATION_COMPLETE:
        altSign = One;
        for (IndexType i=1; i<=n; ++i) {
            x(i) = altSign*(One + T(i-1)/T(n-1));
            altSign = -altSign;
        }
        kase = 1;
        iSave(1) = 5;
        return;
//
//  ................ ENTRY   (ISAVE( 1 ) = 5)
//  X HAS BEEN OVERWRITTEN BY A*X.
//
    case 5:
        const T temp = Two*(blas::asum(x) / T(3*n));
        if (temp>est) {
            v = x;
            est = temp;
        }
    }

QUIT:
    kase = 0;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename  VV, typename VX, typename VSGN, typename EST,
          typename KASE, typename VSAVE>
void
lacn2(DenseVector<VV> &v, DenseVector<VX> &x, DenseVector<VSGN> &sgn,
      EST &est, KASE &kase, DenseVector<VSAVE> &iSave)
{
    typedef typename DenseVector<VV>::IndexType  IndexType;

    cxxlapack::lacn2<IndexType>(v.length(), v.data(), x.data(), sgn.data(),
                                est, kase, iSave.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename  VV, typename VX, typename VSGN, typename EST,
          typename KASE, typename VSAVE>
void
lacn2(DenseVector<VV> &v, DenseVector<VX> &x, DenseVector<VSGN> &sgn,
      EST &est, KASE &kase, DenseVector<VSAVE> &iSave)
{
    typedef typename DenseVector<VV>::IndexType  IndexType;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType n = v.length();

    ASSERT(v.firstIndex()==1);
    ASSERT(v.length()==n);

    ASSERT(x.firstIndex()==1);
    ASSERT(x.length()==n);

    ASSERT(sgn.firstIndex()==1);
    ASSERT(sgn.length()==n);

    ASSERT(iSave.firstIndex()==1);
    ASSERT(iSave.length()==3);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    const typename DenseVector<VV>::NoView      v_org     = v; 
    const typename DenseVector<VX>::NoView      x_org     = x;
    const typename DenseVector<VSGN>::NoView    sgn_org   = sgn;
    const EST                                   est_org   = est;
    const IndexType                             kase_org  = kase; 
    const typename DenseVector<VSAVE>::NoView   iSave_org = iSave;
#   endif

//
//  Call implementation
//
    lacn2_generic(v, x, sgn, est, kase, iSave);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by generic implementation
//
    const typename DenseVector<VV>::NoView      v_generic     = v; 
    const typename DenseVector<VX>::NoView      x_generic     = x;
    const typename DenseVector<VSGN>::NoView    sgn_generic   = sgn;
    const EST                                   est_generic   = est;
    const IndexType                             kase_generic  = kase; 
    const typename DenseVector<VSAVE>::NoView   iSave_generic = iSave;

//
//  restore output arguments
//
    v     = v_org; 
    x     = x_org;
    sgn   = sgn_org;
    est   = est_org;
    kase  = kase_org; 
    iSave = iSave_org;

//
//  Compare generic results with results from the native implementation
//
    external::lacn2(v, x, sgn, est, kase, iSave);

    bool failed = false;
    if (! isIdentical(v_generic, v, "v_generic", "v")) {
        std::cerr << "CXXLAPACK: v_generic = " << v_generic << std::endl;
        std::cerr << "F77LAPACK: v = " << v << std::endl;
        failed = true;
    }
    if (! isIdentical(x_generic, x, "x_generic", "x")) {
        std::cerr << "CXXLAPACK: x_generic = " << x_generic << std::endl;
        std::cerr << "F77LAPACK: x = " << x << std::endl;
        failed = true;
    }
    if (! isIdentical(sgn_generic, sgn, "sgn_generic", "sgn")) {
        std::cerr << "CXXLAPACK: sgn_generic = " << sgn_generic << std::endl;
        std::cerr << "F77LAPACK: sgn = " << sgn << std::endl;
        failed = true;
    }
    if (! isIdentical(est_generic, est, "est_generic", "est")) {
        std::cerr << "CXXLAPACK: est_generic = " << est_generic << std::endl;
        std::cerr << "F77LAPACK: est = " << est << std::endl;
        failed = true;
    }
    if (! isIdentical(kase_generic, kase, "kase_generic", "kase")) {
        std::cerr << "CXXLAPACK: kase_generic = " << kase_generic << std::endl;
        std::cerr << "F77LAPACK: kase = " << kase << std::endl;
        failed = true;
    }
    if (! isIdentical(iSave_generic, iSave, "iSave_generic", "iSave")) {
        std::cerr << "CXXLAPACK: iSave_generic = "
                  << iSave_generic << std::endl;
        std::cerr << "F77LAPACK: iSave = " << iSave << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: lacn2.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: lacn2.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename  VV, typename VX, typename VSGN, typename EST,
          typename KASE, typename VSAVE>
void
lacn2(VV &&v, VX &&x, VSGN &&sgn, EST &&est, KASE &&kase, VSAVE &&iSave)
{
    lacn2(v, x, sgn, est, kase, iSave);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_AUX_LACN2_TCC
