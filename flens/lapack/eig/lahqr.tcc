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
      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_LAHQR_TCC
#define FLENS_LAPACK_EIG_LAHQR_TCC 1

#include <cmath>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr_generic(bool                  wantT,
              bool                  wantZ,
              IndexType             iLo,
              IndexType             iHi,
              GeMatrix<MH>          &H,
              DenseVector<VWR>      &wr,
              DenseVector<VWI>      &wi,
              IndexType             iLoZ,
              IndexType             iHiZ,
              GeMatrix<MZ>          &Z)
{
    using std::abs;
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType  T;

    const Underscore<IndexType>     _;
    const T                         Zero(0), One(1), Two(2);
    const T                         Dat1 = T(3)/T(4),
                                    Dat2 = T(-0.4375);
    const IndexType                 itMax = 30;
    const IndexType                 n = H.numRows();

    typedef typename GeMatrix<MH>::VectorView   VectorView;
    T           vBuffer[3];
    VectorView  v = typename VectorView::Engine(3, vBuffer);

//
//  Quick return if possible
//
    if (n==0) {
        return 0;
    }
    if (iLo==iHi) {
        wr(iLo) = H(iLo, iLo);
        wi(iLo) = Zero;
        return 0;
    }
//
//  ==== clear out the trash ====
//
    for (IndexType j=iLo; j<=iHi-3; ++j) {
        H(j+2, j) = Zero;
        H(j+3, j) = Zero;
    }
    if (iLo<=iHi-2) {
        H(iHi, iHi-2) = Zero;
    }

    const IndexType nh = iHi - iLo + 1;
//
//  Set machine-dependent constants for the stopping criterion.
//
    T safeMin = lamch<T>(SafeMin);
    T safeMax = One / safeMin;
    labad(safeMin, safeMax);
    const T ulp = lamch<T>(Precision);
    const T smallNum = safeMin*(T(nh)/ulp);
//
//  i1 and i2 are the indices of the first row and last column of H
//  to which transformations must be applied. If eigenvalues only are
//  being computed, i1 and i2 are set inside the main loop.
//
    IndexType i1 = -1, i2 = -1;
    if (wantT) {
        i1 = 1;
        i2 = n;
    }
//
//  The main loop begins here. Variable i is the loop index and decreases from
//  iHi to iLo in steps of 1 or 2. Each iteration of the loop works
//  with the active submatrix in rows and columns l to i.
//  Eigenvalues i+1 to iHi have already converged. Either l = iLo or
//  H(l,l-1) is negligible so that the matrix splits.
//
    IndexType i = iHi;
    while (true) {
        IndexType l = iLo;
        if (i<iLo) {
            break;
        }
//
//      Perform QR iterations on rows and columns iLo to i until a
//      submatrix of order 1 or 2 splits off at the bottom because a
//      subdiagonal element has become negligible.
//
        IndexType its;
        for (its=0; its<=itMax; ++its) {
//
//          Look for a single small subdiagonal element.
//
            IndexType k;
            for (k=i; k>=l+1; --k) {
                if (abs(H(k,k-1))<=smallNum) {
                    break;
                }
                T test = abs(H(k-1,k-1)) + abs(H(k,k));
                if (test==Zero) {
                    if (k-2>=iLo) {
                        test += abs(H(k-1,k-2));
                    }
                    if (k+1<=iHi) {
                        test += abs(H(k+1,k));
                    }
                }
//              ==== The following is a conservative small subdiagonal
//              .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
//              .    1997). It has better mathematical foundation and
//              .    improves accuracy in some cases.  ====
                if (abs(H(k,k-1))<=ulp*test) {
                    const T ab = max(abs(H(k, k-1)), abs(H(k-1, k)));
                    const T ba = min(abs(H(k, k-1)), abs(H(k-1, k)));
                    const T aa = max(abs(H(k,   k)), abs(H(k-1, k-1)-H(k, k)));
                    const T bb = min(abs(H(k,   k)), abs(H(k-1, k-1)-H(k, k)));
                    const T s = aa + ab;
                    if (ba*(ab/s)<=max(smallNum, ulp*(bb*(aa/s)))) {
                        break;
                    }
                }
            }
            l = k;

            if (l>iLo) {
//
//              H(l,l-1) is negligible
//
                H(l, l-1) = Zero;
            }
//
//          Exit from loop if a submatrix of order 1 or 2 has split off.
//
            if (l>=i-1) {
                break;
            }
//
//          Now the active submatrix is in rows and columns l to i. If
//          eigenvalues only are being computed, only the active submatrix
//          need be transformed.
//
            if (!wantT) {
                i1 = l;
                i2 = i;
            }

            T H11, H12, H21, H22;
            T rt1r, rt2r, rt1i, rt2i;
            if (its==10) {
//
//              Exceptional shift.
//
                const T s = abs(H(l+1,l)) + abs(H(l+2, l+1));
                H11 = Dat1*s + H(l,l);
                H12 = Dat2*s;
                H21 = s;
                H22 = H11;
            } else if (its==20) {
//
//              Exceptional shift.
//
                const T s = abs(H(i,i-1)) + abs(H(i-1,i-2));
                H11 = Dat1*s + H(i,i);
                H12 = Dat2*s;
                H21 = s;
                H22 = H11;
            } else {
//
//              Prepare to use Francis' double shift
//              (i.e. 2nd degree generalized Rayleigh quotient)
//
                H11 = H(i-1, i-1);
                H21 = H(  i, i-1);
                H12 = H(i-1,   i);
                H22 = H(  i,   i);
            }

            const T s = abs(H11) + abs(H12) + abs(H21) + abs(H22);
            if (s==Zero) {
                rt1r = Zero;
                rt1i = Zero;
                rt2r = Zero;
                rt2i = Zero;
            } else {
                H11 /= s;
                H21 /= s;
                H12 /= s;
                H22 /= s;
                const T tr = (H11+H22) / Two;
                const T det = (H11-tr)*(H22-tr) - H12*H21;
                const T rtDisc = sqrt(abs(det));
                if (det>=Zero) {
//
//                  ==== complex conjugate shifts ====
//
                    rt1r = tr*s;
                    rt2r = rt1r;
                    rt1i = rtDisc*s;
                    rt2i = -rt1i;
                } else {
//
//                  ==== real shifts (use only one of them)  ====
//
                    rt1r = tr + rtDisc;
                    rt2r = tr - rtDisc;
                    if (abs(rt1r-H22)<=abs(rt2r-H22)) {
                        rt1r *= s;
                        rt2r = rt1r;
                    } else {
                        rt2r *= s;
                        rt1r = rt2r;
                    }
                    rt1i = Zero;
                    rt2i = Zero;
                }
            }
//
//          Look for two consecutive small subdiagonal elements.
//
            IndexType m;
            for (m=i-2; m>=l; --m) {
//              Determine the effect of starting the double-shift QR
//              iteration at row m, and see if this would make H(m,m-1)
//              negligible.  (The following uses scaling to avoid
//              overflows and most underflows.)
//
                T H21S = H(m+1,m);
                T s = abs(H(m,m)-rt2r) + abs(rt2i) + abs(H21S);

                H21S = H(m+1,m)/s;
                v(1) = H21S*H(m,m+1)
                     + (H(m,m)-rt1r)*((H(m,m)-rt2r)/s)
                     - rt1i*(rt2i/s);
                v(2) = H21S*(H(m,m) + H(m+1,m+1)-rt1r-rt2r);
                v(3) = H21S*H(m+2,m+1);

                s = abs(v(1)) + abs(v(2)) + abs(v(3));
                v(1) /= s;
                v(2) /= s;
                v(3) /= s;

                if (m==l) {
                    break;
                }

                const T value1 = abs(H(m,m-1))*(abs(v(2))+abs(v(3)));
                const T value2 = ulp*abs(v(1))
                                 *(abs(H(m-1,m-1))+abs(H(m,m))+abs(H(m+1,m+1)));
                if (value1<=value2) {
                    break;
                }
            }
//
//          Double-shift QR step
//
            for (k=m; k<=i-1; ++k) {
//
//              The first iteration of this loop determines a reflection G
//              from the vector v and applies it from left and right to H,
//              thus creating a nonzero bulge below the subdiagonal.
//
//              Each subsequent iteration determines a reflection G to
//              restore the Hessenberg form in the (k-1)th column, and thus
//              chases the bulge one step toward the bottom of the active
//              submatrix. 'nr' is the order of G.
//
                T   t1, t2, t3, v2, v3;

                IndexType nr = min(IndexType(3), i-k+1);
                if (k>m) {
                    v(_(1,nr)) = H(_(k,k+nr-1),k-1);
                }
                larfg(nr, v(1), v(_(2,nr)), t1);
                if (k>m) {
                    H(k,  k-1) = v(1);
                    H(k+1,k-1) = Zero;
                    if (k<i-1) {
                        H(k+2,k-1) = Zero;
                    }
                } else if (m>l) {
//                  ==== Use the following instead of
//                  .    H( K, K-1 ) = -H( K, K-1 ) to
//                  .    avoid a bug when v(2) and v(3)
//                  .    underflow. ====
                    H(k, k-1) *= (One-t1);
                }
                v2 = v(2);
                t2 = t1*v2;
                if (nr==3) {
                    v3 = v(3);
                    t3 = t1*v3;
//
//                  Apply G from the left to transform the rows of the matrix
//                  in columns k to i2.
//
                    for (IndexType j=k; j<=i2; ++j) {
                        const T sum = H(k,j) + v2*H(k+1,j) + v3*H(k+2,j);
                        H(k,   j) -= sum*t1;
                        H(k+1, j) -= sum*t2;
                        H(k+2, j) -= sum*t3;
                    }
//
//                  Apply G from the right to transform the columns of the
//                  matrix in rows i1 to min(k+3,i).
//
                    for (IndexType j=i1; j<=min(k+3,i); ++j) {
                        const T sum = H(j, k) + v2*H(j,k+1) + v3*H(j,k+2);
                        H(j, k  ) -= sum*t1;
                        H(j, k+1) -= sum*t2;
                        H(j, k+2) -= sum*t3;
                    }

                    if (wantZ) {
//
//                      Accumulate transformations in the matrix Z
//
                        for (IndexType j=iLoZ; j<=iHiZ; ++j) {
                            const T sum = Z(j, k) + v2*Z(j, k+1) + v3*Z(j, k+2);
                            Z(j, k  ) -= sum*t1;
                            Z(j, k+1) -= sum*t2;
                            Z(j, k+2) -= sum*t3;
                        }
                    }
                } else if (nr==2) {
//
//                  Apply G from the left to transform the rows of the matrix
//                  in columns K to I2.
//
                    for (IndexType j=k; j<=i2; ++j) {
                        const T sum = H(k, j) + v2*H(k+1, j);
                        H(k,   j) -= sum*t1;
                        H(k+1, j) -= sum*t2;
                    }
//
//                  Apply G from the right to transform the columns of the
//                  matrix in rows i1 to min(k+3,i).
//
                    for (IndexType j=i1; j<=i; ++j) {
                        const T sum = H(j, k) + v2*H(j, k+1);
                        H(j, k  ) -= sum*t1;
                        H(j, k+1) -= sum*t2;
                    }

                    if (wantZ) {
//
//                      Accumulate transformations in the matrix Z
//
                        for (IndexType j=iLoZ; j<=iHiZ; ++j) {
                            const T sum = Z(j, k) + v2*Z(j, k+1);
                            Z(j, k  ) -= sum*t1;
                            Z(j, k+1) -= sum*t2;
                        }
                    }
                }
            }
        }
//
//      Failure to converge in remaining number of iterations
//
        if (its>itMax) {
            return i;
        }

        if (l==i) {
//
//          H(I,I-1) is negligible: one eigenvalue has converged.
//
            wr(i) = H(i, i);
            wi(i) = Zero;
        } else if (l==i-1) {
//
//          H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
//
//          Transform the 2-by-2 submatrix to standard Schur form,
//          and compute and store the eigenvalues.
//
            T cs, sn;
            lanv2(H(i-1,i-1), H(i-1,i), H(i,i-1), H(i,i),
                  wr(i-1), wi(i-1), wr(i), wi(i),
                  cs, sn);

            if (wantT) {
//
//              Apply the transformation to the rest of H.
//
                if (i2>i) {
                    const auto cols = _(i+1,i2);
                    blas::rot(H(i-1,cols), H(i,cols), cs, sn);
                }
                const auto rows = _(i1,i-2);
                blas::rot(H(rows, i-1), H(rows, i), cs, sn);
            }
            if (wantZ) {
//
//              Apply the transformation to Z.
//
                const auto rows = _(iLoZ, iHiZ);
                blas::rot(Z(rows, i-1), Z(rows, i), cs, sn);
            }
        }
//
//      return to start of the main loop with new value of I.
//
        i = l - 1;
    }
    return 0;
}

//== interface for native lapack ===============================================

#ifdef TODO_CHECK_CXXLAPACK

template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr_native(bool                   wantT,
             bool                   wantZ,
             IndexType              iLo,
             IndexType              iHi,
             GeMatrix<MH>           &H,
             DenseVector<VWR>       &wr,
             DenseVector<VWI>       &wi,
             IndexType              iLoZ,
             IndexType              iHiZ,
             GeMatrix<MZ>           &Z)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    const LOGICAL    WANTT  = wantT;
    const LOGICAL    WANTZ  = wantZ;
    const INTEGER    N      = H.numRows();
    const INTEGER    ILO    = iLo;
    const INTEGER    IHI    = iHi;
    const INTEGER    LDH    = H.leadingDimension();
    const INTEGER    ILOZ   = iLoZ;
    const INTEGER    IHIZ   = iHiZ;
    const INTEGER    LDZ    = Z.leadingDimension();
    INTEGER          INFO;

    if (IsSame<T,DOUBLE>::value) {
        LAPACK_IMPL(dlahqr)(&WANTT,
                            &WANTZ,
                            &N,
                            &ILO,
                            &IHI,
                            H.data(),
                            &LDH,
                            wr.data(),
                            wi.data(),
                            &ILOZ,
                            &IHIZ,
                            Z.data(),
                            &LDZ,
                            &INFO);
    } else {
        ASSERT(0);
    }
    return INFO;
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr(bool                  wantT,
      bool                  wantZ,
      IndexType             iLo,
      IndexType             iHi,
      GeMatrix<MH>          &H,
      DenseVector<VWR>      &wr,
      DenseVector<VWI>      &wi,
      IndexType             iLoZ,
      IndexType             iHiZ,
      GeMatrix<MZ>          &Z)
{
    LAPACK_DEBUG_OUT("lahqr");

//
//  Test the input parameters
//
    using std::max;

    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());
    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==H.numRows());
    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==H.numRows());
    ASSERT(wantZ || (Z.numRows()==H.numCols()));
    ASSERT(wantZ || (Z.numCols()==H.numCols()));

    // 1 <= ILO <= max(1,IHI); IHI <= N.
    ASSERT(1<=iLo);
    ASSERT(iLo<=max(IndexType(1), iHi));
    ASSERT(iHi<=H.numRows());

    // 1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
    ASSERT(1<=iLoZ);
    ASSERT(iLoZ<=iLo);
    ASSERT(iHi<=iHiZ);
    ASSERT(iHiZ<=H.numRows());

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MH>::NoView          H_org    = H;
    typename GeMatrix<MH>::NoView          Z_org    = Z;

    typename GeMatrix<MH>::NoView          _H       = H;
    typename DenseVector<VWR>::NoView      _wr      = wr;
    typename DenseVector<VWI>::NoView      _wi      = wi;
    typename GeMatrix<MZ>::NoView          _Z       = Z;
#   endif

//
//  Call implementation
//
    IndexType info = lahqr_generic(wantT, wantZ, iLo,  iHi,
                                   H, wr, wi, iLoZ, iHiZ, Z);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    IndexType _info = lahqr_native(wantT, wantZ, iLo,  iHi,
                                   _H, _wr, _wi, iLoZ, iHiZ, _Z);

    bool failed = false;
    if (! isIdentical(H, _H, " H", "_H")) {
        std::cerr << "CXXLAPACK:  H = " << H << std::endl;
        std::cerr << "F77LAPACK: _H = " << _H << std::endl;
        failed = true;
    }

    if (! isIdentical(wr, _wr, " wr", "_wr")) {
        std::cerr << "CXXLAPACK:  wr = " << wr << std::endl;
        std::cerr << "F77LAPACK: _wr = " << _wr << std::endl;
        failed = true;
    }

    if (! isIdentical(wi, _wi, " wi", "_wi")) {
        std::cerr << "CXXLAPACK:  wi = " << wi << std::endl;
        std::cerr << "F77LAPACK: _wi = " << _wi << std::endl;
        failed = true;
    }

    if (! isIdentical(Z, _Z, " Z", "_Z")) {
        std::cerr << "CXXLAPACK:  Z = " << Z << std::endl;
        std::cerr << "F77LAPACK: _Z = " << _Z << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "H_org = " << H_org << std::endl;
        std::cerr << "Z_org = " << Z_org << std::endl;
        std::cerr << "wantT = " << wantT
                  << ", wantZ = " << wantZ
                  << ", iLo = " << iLo
                  << ", iHi = " << iHi
                  << std::endl;
        std::cerr << "error in: lahqr.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: lahqr.tcc" << std::endl;
    }
#   endif
    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr(bool          wantT,
      bool          wantZ,
      IndexType     iLo,
      IndexType     iHi,
      MH            &&H,
      VWR           &&wr,
      VWI           &&wi,
      IndexType     iLoZ,
      IndexType     iHiZ,
      MZ            &&Z)
{
    return lahqr(wantT, wantZ, iLo, iHi, H, wr, wi, iLoZ, iHiZ, Z);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LAHQR_TCC
