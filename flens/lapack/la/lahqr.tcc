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
      SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAHQR_TCC
#define FLENS_LAPACK_LA_LAHQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr_impl(bool                  wantT,
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

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ>
IndexType
lahqr_impl(bool                  wantT,
           bool                  wantZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VW>       &w,
           IndexType             iLoZ,
           IndexType             iHiZ,
           GeMatrix<MZ>          &Z)
{
    using cxxblas::pow;
    using std::abs;
    using std::conj;
    using std::imag;
    using std::max;
    using std::min;
    using std::real;
    using cxxblas::abs1;

    typedef typename GeMatrix<MH>::ElementType          T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;

    const Underscore<IndexType>     _;
    const T                         Zero(0), One(1);
    const PT                        RZero(0), ROne(1), RHalf(0.5);

    const PT                        Dat1 = PT(3)/PT(4);
    const IndexType                 itMax = 30;
    const IndexType                 n = H.numRows();


    typedef typename GeMatrix<MH>::VectorView           VectorView;

    T           vBuffer[2];
    VectorView  v = typename VectorView::Engine(2, vBuffer);

//
//  Quick return if possible
//
    if (n==0) {
        return 0;
    }
    if (iLo==iHi) {
        w(iLo) = H(iLo, iLo);
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

//
//  ==== ensure that subdiagonal entries are real ====
//
    const IndexType jLo = (wantT) ? 1 : iLo;
    const IndexType jHi = (wantT) ? n : iHi;

    for (IndexType i=iLo+1; i<=iHi; ++i) {
        if (imag(H(i,i-1))!=RZero) {
//          ==== The following redundant normalization
//          .    avoids problems with both gradual and
//          .    sudden underflow in ABS(H(I,I-1)) ====
            T scale = H(i,i-1) / abs1(H(i,i-1));
            scale = conj(scale) / abs(scale);
            H(i,i-1) = abs(H(i,i-1));
            H(i, _(i,jHi))            *= scale;
            H(_(jLo,min(jHi,i+1)), i) *= conj(scale);

            if (wantZ) {
                Z(_(iLoZ,iHiZ),i) *= conj(scale);
            }
        }
    }

    const IndexType nh = iHi - iLo + 1;
//
//  Set machine-dependent constants for the stopping criterion.
//
    PT safeMin = lamch<PT>(SafeMin);
    PT safeMax = ROne / safeMin;
    labad(safeMin, safeMax);
    const PT ulp      = lamch<PT>(Precision);
    const PT smallNum = safeMin*(PT(nh)/ulp);
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
        if (i<iLo) {
            break;
        }
//
//      Perform QR iterations on rows and columns iLo to i until a
//      submatrix of order 1 or 2 splits off at the bottom because a
//      subdiagonal element has become negligible.
//
        IndexType l = iLo;
        IndexType its;
        for (its=0; its<=itMax; ++its) {
//
//          Look for a single small subdiagonal element.
//
            IndexType k;
            for (k=i; k>=l+1; --k) {
                if (abs1(H(k,k-1))<=smallNum) {
                    break;
                }
                PT test = abs1(H(k-1,k-1)) + abs1(H(k,k));
                if (test==RZero) {
                    if (k-2>=iLo) {
                        test += abs(real(H(k-1,k-2)));
                    }
                    if (k+1<=iHi) {
                        test += abs(real(H(k+1,k)));
                    }
                }
//              ==== The following is a conservative small subdiagonal
//              .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
//              .    1997). It has better mathematical foundation and
//              .    improves accuracy in some cases.  ====
                if (abs(real(H(k,k-1)))<=ulp*test) {
                    const PT ab = max(abs1(H(k, k-1)),
                                      abs1(H(k-1, k)));
                    const PT ba = min(abs1(H(k, k-1)),
                                      abs1(H(k-1, k)));
                    const PT aa = max(abs1(H(k, k)),
                                      abs1(H(k-1, k-1)-H(k, k)));
                    const PT bb = min(abs1(H(k, k)),
                                      abs1(H(k-1, k-1)-H(k, k)));
                    const PT s = aa + ab;
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
            if (l>=i) {
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

            T   H11, H11S, H22, t, u, x, y;
            PT  H10, H21, s, sx;
            if (its==10) {
//
//              Exceptional shift.
//
                s = Dat1*abs(real(H(l+1,l)));
                t = s + H(l,l);
            } else if (its==20) {
//
//              Exceptional shift.
//
                s = Dat1*abs(real(H(i,i-1)));
                t = s + H(i,i);
            } else {
//
//
//              Wilkinson's shift.
//
//
                t = H(i,i);
                u = sqrt(H(i-1,i)) * sqrt(H(i,i-1));
                s = abs1(u);

                if (s != RZero) {
                    x  = RHalf * (H(i-1,i-1) - t);
                    sx = abs1(x);
                    s  = max(s, abs1(x));
                    y  = s * sqrt(cxxblas::pow(x/s,2) + cxxblas::pow(u/s,2));
                    if (sx > RZero) {
                        PT tmp = real(x/sx)*real(y)
                               + imag(x/sx)*imag(y);

                        if (tmp<RZero) {
                            y = -y;
                        }
                    }
                    t -= u*ladiv(u, x+y);
                }
            }
//
//          Look for two consecutive small subdiagonal elements.
//
            IndexType m;
            for (m=i-1; m>=l+1; --m) {
//
//          Determine the effect of starting the single-shift QR
//          iteration at row M, and see if this would make H(M,M-1)
//          negligible.
//
                H11 = H(m, m);
                H22 = H(m+1, m+1);
                H11S = H11 - t;
                H21 = real(H(m+1,m));
                s = abs1(H11S) + abs(H21);
                H11S /= s;
                H21  /= s;
                v(1) = H11S;
                v(2) = H21;
                H10 = real(H(m,m-1));
                if (abs(H10)*abs(H21)<=ulp*(abs1(H11S)*(abs1(H11)+abs1(H22)))) {
                    goto SINGLE_SHIFT;
                }
            }
            H11  = H(l, l);
            H22  = H(l+1, l+1);
            H11S = H11 - t;
            H21  = real(H(l+1,l));
            s    = abs1(H11S) + abs(H21);
            H11S /= s;
            H21  /= s;
            v(1) = H11S;
            v(2) = H21;

        SINGLE_SHIFT:
//
//          Single-shift QR step
//
            for (k=m; k<=i-1; ++k) {
//
//              The first iteration of this loop determines a reflection G
//              from the vector V and applies it from left and right to H,
//              thus creating a nonzero bulge below the subdiagonal.
//
//              Each subsequent iteration determines a reflection G to
//              restore the Hessenberg form in the (K-1)th column, and thus
//              chases the bulge one step toward the bottom of the active
//              submatrix.
//
//              V(2) is always real before the call to ZLARFG, and hence
//              after the call T2 ( = T1*V(2) ) is also real.
//
                T  t1, v2;
                PT t2;

                if (k>m) {
                    v = H(_(k,k+1),k-1);
                }
                larfg(2, v(1), v(_(2,2)), t1);
                if (k>m) {
                    H(k,  k-1) = v(1);
                    H(k+1,k-1) = Zero;
                }
                v2 = v(2);
                t2 = real(t1*v2);
//
//              Apply G from the left to transform the rows of the matrix
//              in columns K to I2.
//
                for (IndexType j=k; j<=i2; ++j) {
                    const T sum = conj(t1)*H(k,j) + t2*H(k+1,j);
                    H(k,  j) -= sum;
                    H(k+1,j) -= sum*v2;
                }
//
//              Apply G from the right to transform the columns of the
//              matrix in rows I1 to min(K+2,I).
//
                for (IndexType j=i1; j<=min(k+2,i); ++j) {
                    T sum = t1*H(j,k) + t2*H(j,k+1);
                    H(j,k)   -= sum;
                    H(j,k+1) -= sum*conj(v2);
                }

                if (wantZ) {
//
//                  Accumulate transformations in the matrix Z
//
                    for (IndexType j=iLoZ; j<=iHiZ; ++j) {
                        T sum = t1*Z(j,k) + t2*Z(j,k+1);
                        Z(j,k)   -= sum;
                        Z(j,k+1) -= sum*conj(v2);
                    }
                }

                if (k==m && m>l) {
//
//                  If the QR step was started at row M > L because two
//                  consecutive small subdiagonals were found, then extra
//                  scaling must be performed to ensure that H(M,M-1) remains
//                  real.
//
                    T tmp = One - t1;
                    tmp /= abs(tmp);
                    H(m+1,m) *= conj(tmp);
                    if (m+2 <= i) {
                        H(m+2,m+1) *= tmp;
                    }
                    for (IndexType j=m; j<=i; ++j) {
                        if (j != m+1) {
                            if (i2>j) {
                                H(j,_(j+1,i2)) *= tmp;
                            }
                            H(_(i1,j-1),j) *= conj(tmp);
                            if (wantZ) {
                                Z(_(iLoZ,iHiZ),j) *= conj(tmp);
                            }
                        }
                    }
                }
            }

//
//          Ensure that H(I,I-1) is real.
//
            T tmp = H(i,i-1);
            if (imag(tmp) != RZero) {
                PT rTmp = abs(tmp);
                H(i,i-1) = rTmp;
                tmp /= rTmp;
                if (i2>i) {
                    H(i,_(i+1,i2)) *= conj(tmp);
                }
                H(_(i1,i-1),i) *= tmp;
                if (wantZ) {
                    Z(_(iLoZ,iHiZ),i) *= tmp;
                }
            }
        }
//
//      Failure to converge in remaining number of iterations
//
        if (its>itMax) {
            return i;
        }

//
//      H(I,I-1) is negligible: one eigenvalue has converged.
//
        w(i) = H(i,i);
//
//      return to start of the main loop with new value of I.
//
        i = l - 1;
    }
    return 0;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
IndexType
lahqr_impl(bool                   wantT,
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
    IndexType  info;
    info = cxxlapack::lahqr<IndexType>(wantT,
                                       wantZ,
                                       H.numRows(),
                                       iLo,
                                       iHi,
                                       H.data(),
                                       H.leadingDimension(),
                                       wr.data(),
                                       wi.data(),
                                       iLoZ,
                                       iHiZ,
                                       Z.data(),
                                       Z.leadingDimension());
    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ>
IndexType
lahqr_impl(bool                   wantT,
           bool                   wantZ,
           IndexType              iLo,
           IndexType              iHi,
           GeMatrix<MH>           &H,
           DenseVector<VW>        &w,
           IndexType              iLoZ,
           IndexType              iHiZ,
           GeMatrix<MZ>           &Z)
{
    IndexType  info;
    info = cxxlapack::lahqr<IndexType>(wantT,
                                       wantZ,
                                       H.numRows(),
                                       iLo,
                                       iHi,
                                       H.data(),
                                       H.leadingDimension(),
                                       w.data(),
                                       iLoZ,
                                       iHiZ,
                                       Z.data(),
                                       Z.leadingDimension());
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ>
typename RestrictTo<IsRealGeMatrix<MH>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MZ>::value,
         IndexType>::Type
lahqr(bool              wantT,
      bool              wantZ,
      IndexType         iLo,
      IndexType         iHi,
      MH                &&H,
      VWR               &&wr,
      VWI               &&wi,
      IndexType         iLoZ,
      IndexType         iHiZ,
      MZ                &&Z)
{
    LAPACK_DEBUG_OUT("lahqr");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type        MatrixH;
    typedef typename RemoveRef<VWR>::Type       VectorWR;
    typedef typename RemoveRef<VWI>::Type       VectorWI;
    typedef typename RemoveRef<MZ>::Type        MatrixZ;
#   endif

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
    typename MatrixH::NoView        H_org    = H;
    typename MatrixZ::NoView        Z_org    = Z;

    typename MatrixH::NoView        _H       = H;
    typename VectorWR::NoView       _wr      = wr;
    typename VectorWI::NoView       _wi      = wi;
    typename MatrixZ::NoView        _Z       = Z;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::lahqr_impl(wantT, wantZ, iLo,  iHi,
                                               H, wr, wi, iLoZ, iHiZ, Z);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    IndexType _info = external::lahqr_impl(wantT, wantZ, iLo,  iHi,
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

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ>
typename RestrictTo<IsComplexGeMatrix<MH>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MZ>::value,
         IndexType>::Type
lahqr(bool              wantT,
      bool              wantZ,
      IndexType         iLo,
      IndexType         iHi,
      MH                &&H,
      VW                &&w,
      IndexType         iLoZ,
      IndexType         iHiZ,
      MZ                &&Z)
{
    LAPACK_DEBUG_OUT("lahqr (complex)");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type        MatrixH;
    typedef typename RemoveRef<VW>::Type        VectorW;
    typedef typename RemoveRef<MZ>::Type        MatrixZ;
#   endif

//
//  Test the input parameters
//
    using std::max;

    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());
    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==H.numRows());
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
    typename MatrixH::NoView        H_org    = H;
    typename MatrixZ::NoView        Z_org    = Z;

    typename MatrixH::NoView        _H       = H;
    typename VectorW::NoView        _w       = w;
    typename MatrixZ::NoView        _Z       = Z;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::lahqr_impl(wantT, wantZ, iLo,  iHi,
                                               H, w, iLoZ, iHiZ, Z);

//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    IndexType _info = external::lahqr_impl(wantT, wantZ, iLo,  iHi,
                                           _H, _w, iLoZ, iHiZ, _Z);

    bool failed = false;
    if (! isIdentical(H, _H, " H", "_H")) {
        std::cerr << "CXXLAPACK:  H = " << H << std::endl;
        std::cerr << "F77LAPACK: _H = " << _H << std::endl;
        failed = true;
    }

    if (! isIdentical(w, _w, " w", "_w")) {
        std::cerr << "CXXLAPACK:  w = " << w << std::endl;
        std::cerr << "F77LAPACK: _w = " << _w << std::endl;
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

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAHQR_TCC
