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
      SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
     $                   LDT, NV, WV, LDWV, WORK, LWORK )
      SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
     $                   NV, WV, LDWV, WORK, LWORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2.2)                        --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *  -- June 2010                                                       --
 *
 */

#ifndef FLENS_LAPACK_LA_LAQR2_TCC
#define FLENS_LAPACK_LA_LAQR2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Workspace query (real variant)
//
template <typename IndexType, typename MT>
typename RestrictTo<IsReal<typename MT::ElementType>::value,
         IndexType>::Type
laqr2_wsq_impl(IndexType                 kTop,
               IndexType                 kBot,
               IndexType                 nw,
               const GeMatrix<MT>        &T)
{
    using std::max;
    using std::min;

    const Underscore<IndexType> _;

//
//  ==== Estimate optimal workspace. ====
//
    IndexType jw = min(nw, kBot-kTop+1);
    auto  _T  = T(_(1,jw),_(1,jw));

    IndexType lWorkOpt;
    if (jw<=2) {
        lWorkOpt = 1;
    } else {
//
//      ==== Workspace query call to DGEHRD ====
//
        IndexType lWork1 = hrd_wsq(IndexType(1), jw-1, _T);
//
//      ==== Workspace query call to DORMHR ====
//
        IndexType lWork2 = ormhr_wsq(Right, NoTrans, IndexType(1), jw-1, _T);
//
//      ==== Optimal workspace ====
//
        lWorkOpt = jw+max(lWork1,lWork2);
    }
    return lWorkOpt;
}

//
//  Workspace query (complex variant)
//
template <typename IndexType, typename MT>
typename RestrictTo<IsComplex<typename MT::ElementType>::value,
         IndexType>::Type
laqr2_wsq_impl(IndexType                 kTop,
               IndexType                 kBot,
               IndexType                 nw,
               const GeMatrix<MT>        &T)
{
    using std::max;
    using std::min;

    const Underscore<IndexType> _;

//
//  ==== Estimate optimal workspace. ====
//
    IndexType jw = min(nw, kBot-kTop+1);
    auto  _T  = T(_(1,jw),_(1,jw));

    IndexType lWorkOpt;
    if (jw<=2) {
        lWorkOpt = 1;
    } else {
//
//      ==== Workspace query call to DGEHRD ====
//
        IndexType lWork1 = hrd_wsq(IndexType(1), jw-1, _T);
//
//      ==== Workspace query call to DORMHR ====
//
        IndexType lWork2 = unmhr_wsq(Right, NoTrans, IndexType(1), jw-1, _T);
//
//      ==== Optimal workspace ====
//
        lWorkOpt = jw+max(lWork1,lWork2);
    }
    return lWorkOpt;
}

//
//  Real variant
//
template <typename IndexType, typename MH, typename MZ, typename VSR,
          typename VSI, typename MV, typename MT, typename MWV, typename VWORK>
void
laqr2_impl(bool                      wantT,
           bool                      wantZ,
           IndexType                 kTop,
           IndexType                 kBot,
           IndexType                 nw,
           GeMatrix<MH>              &H,
           IndexType                 iLoZ,
           IndexType                 iHiZ,
           GeMatrix<MZ>              &Z,
           IndexType                 &ns,
           IndexType                 &nd,
           DenseVector<VSR>          &sr,
           DenseVector<VSI>          &si,
           GeMatrix<MV>              &V,
           GeMatrix<MT>              &T,
           GeMatrix<MWV>             &WV,
           DenseVector<VWORK>        &work)
{
    using std::abs;
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType  ElementType;

    const Underscore<IndexType> _;

    const ElementType   Zero(0), One(1);
    const IndexType     n   = H.numRows();
    const IndexType     nh  = T.numCols();
    const IndexType     nv  = WV.numRows();
//
//  ==== Estimate optimal workspace. ====
//
    IndexType   jw      = min(nw, kBot-kTop+1);
    auto        T_jw    = T(_(1,jw),_(1,jw));
    auto        V_jw    = V(_(1,jw),_(1,jw));

    IndexType lWorkOpt;
    if (jw<=2) {
        lWorkOpt = 1;
    } else {
//
//      ==== Workspace query call to DGEHRD ====
//
        IndexType lWork1 = hrd_wsq(IndexType(1), jw-1, T_jw);
//
//      ==== Workspace query call to DORMHR ====
//
        IndexType lWork2 = ormhr_wsq(Right, NoTrans, IndexType(1), jw-1, V_jw);
//
//      ==== Optimal workspace ====
//
        lWorkOpt = jw+max(lWork1,lWork2);
    }
//
//  ==== Apply worksize query ====
//
    if (work.length()==0) {
        work.resize(lWorkOpt);
    }
    const IndexType lWork = work.length();
//
//  ==== Nothing to do ...
//  ... for an empty active block ... ====
    ns = 0;
    nd = 0;
    work(1) = One;
    if (kTop>kBot) {
        return;
    }
//   ... nor for an empty deflation window. ====
    if (nw<1) {
        return;
    }
//
//  ==== Machine constants ====
//
    ElementType safeMin   = lamch<ElementType>(SafeMin);
    ElementType safeMax   = One / safeMin;
    labad(safeMin, safeMax);
    const ElementType ulp      = lamch<ElementType>(Precision);
    const ElementType smallNum = safeMin*(ElementType(n)/ulp);
//
//  ==== Setup deflation window ====
//
    IndexType   kwTop   = kBot - jw + 1;
    auto        sr_kw   = sr(_(kwTop,kBot));
    auto        si_kw   = si(_(kwTop,kBot));
    auto        H_kw    = H(_(kwTop,kBot),_(kwTop,kBot));
    ElementType s;

    if( kwTop==kTop) {
        s = Zero;
    } else {
        s = H(kwTop, kwTop-1);
    }

    if (kBot==kwTop) {
//
//      ==== 1-by-1 deflation window: not much to do ====
//
        sr(kwTop) = H(kwTop,kwTop);
        si(kwTop) = Zero;
        ns = 1;
        nd = 0;
        if (abs(s)<=max(smallNum,ulp*abs(H(kwTop,kwTop)))) {
            ns = 0;
            nd = 1;
            if (kwTop>kTop) {
                H(kwTop, kwTop-1) = Zero;
            }
        }
        work(1) = One;
        return;
    }
//
//  ==== Convert to spike-triangular form.  (In case of a
//  .    rare QR failure, this routine continues to do
//  .    aggressive early deflation using that part of
//  .    the deflation window that converged using INFQR
//  .    here and there to keep track.) ====
//
    T_jw.upper()  = H_kw.upper();
    T_jw.diag(-1) = H_kw.diag(-1);

    V_jw          = Zero;
    V_jw.diag(0)  = One;

    IndexType infoQR;
    infoQR = lahqr(true, true,
                   IndexType(1), jw, T_jw,
                   sr_kw, si_kw,
                   IndexType(1), jw, V_jw);
//
//  ==== DTREXC needs a clean margin near the diagonal ====
//
    for (IndexType j=1; j<=jw-3; ++j) {
        T(j+2, j) = Zero;
        T(j+3, j) = Zero;
    }
    if (jw>2) {
        T(jw,jw-2) = Zero;
    }
//
//  ==== Deflation detection loop ====
//
    ns = jw;
    IndexType iFirst;
    IndexType iLast = infoQR + 1;
    while (iLast<=ns) {
        bool bulge;
        if (ns==1) {
            bulge = false;
        } else {
            bulge = T(ns,ns-1)!=Zero;
        }
//
//      ==== Small spike tip test for deflation ====
//
        if (!bulge) {
//
//          ==== Real eigenvalue ====
//
            ElementType foo = abs(T(ns,ns));
            if (foo==Zero) {
                foo = abs(s);
            }
            if (abs(s*V(1,ns))<=max(smallNum,ulp*foo)) {
//
//               ==== Deflatable ====
//
                ns = ns - 1;
            } else {
//
//               ==== Undeflatable.   Move it up out of the way.
//               .    (DTREXC can not fail in this case.) ====
//
                iFirst = ns;
                trexc(true, T_jw, V_jw, iFirst, iLast, work(_(1,jw)));
                ++iLast;
            }
        } else {
//
//          ==== Complex conjugate pair ====
//
            ElementType foo = abs(T(ns,ns))
                            + sqrt(abs(T(ns,ns-1)))*sqrt(abs(T(ns-1,ns)));
            if (foo==Zero) {
                foo = abs(s);
            }
            if (max(abs(s*V(1,ns)), abs(s*V(1,ns-1)))<=max(smallNum,ulp*foo)) {
//
//              ==== Deflatable ====
//
                ns -= 2;
            } else {
//
//              ==== Undeflatable. Move them up out of the way.
//              .    Fortunately, DTREXC does the right thing with
//              .    ILST in case of a rare exchange failure. ====
//
                iFirst = ns;
                trexc(true, T_jw, V_jw, iFirst, iLast, work(_(1,jw)));
                iLast += 2;
            }
        }
//
//      ==== End deflation detection loop ====
//
    }
//
//     ==== Return to Hessenberg form ====
//
    if (ns==0) {
        s = Zero;
    }

    if (ns<jw) {
//
//     ==== sorting diagonal blocks of T improves accuracy for
//     .    graded matrices.  Bubble sort deals well with
//     .    exchange failures. ====
//
        bool sorted = false;
        IndexType i = ns + 1;
        while (!sorted) {
            sorted = true;

            IndexType kEnd = i - 1;
            i = infoQR + 1;

            IndexType k;
            if (i==ns) {
                k = i + 1;
            } else if (T(i+1,i)==Zero) {
                k = i + 1;
            } else {
                k = i + 2;
            }
            while (k<=kEnd) {
                ElementType evi, evk;

                if (k==i+1) {
                    evi = abs(T(i,i));
                } else {
                    evi = abs(T(i,i)) + sqrt(abs(T(i+1,i)))*sqrt(abs(T(i,i+1)));
                }

                if (k==kEnd) {
                    evk = abs(T(k,k));
                } else if(T(k+1,k)==Zero) {
                    evk = abs(T(k,k));
                } else {
                    evk = abs(T(k,k)) + sqrt(abs(T(k+1,k)))*sqrt(abs(T(k,k+1)));
                }

                if (evi>=evk) {
                    i = k;
                } else {
                    sorted = false;
                    iFirst = i;
                    iLast = k;
                    IndexType info = trexc(true, T_jw, V_jw, iFirst, iLast,
                                           work(_(1,jw)));
                    if (info==0) {
                        i = iLast;
                    } else {
                        i = k;
                    }
                }
                if (i==kEnd) {
                    k = i + 1;
                } else if (T(i+1,i)==Zero) {
                    k = i + 1;
                } else {
                    k = i + 2;
                }
            }
        }
    }
//
//  ==== Restore shift/eigenvalue array from T ====
//
    IndexType i = jw;
    while (i>=infoQR+1) {
        if (i==infoQR+1) {
            sr(kwTop+i-1) = T(i,i);
            si(kwTop+i-1) = Zero;
            i = i - 1;
        } else if (T(i,i-1)==Zero) {
            sr(kwTop+i-1 ) = T(i,i);
            si(kwTop+i-1 ) = Zero;
            i = i - 1;
        } else {
            ElementType aa = T(i-1,i-1);
            ElementType cc = T(i,  i-1);
            ElementType bb = T(i-1,i  );
            ElementType dd = T(i,  i  );
            ElementType cs, sn;
            lanv2(aa, bb, cc, dd,
                  sr(kwTop+i-2), si(kwTop+i-2),
                  sr(kwTop+i-1), si(kwTop+i-1),
                  cs, sn);
            i -= 2;
        }
    }

    if (ns<jw || s==Zero) {

        if (ns>1 && s!=Zero) {
//
//          ==== Reflect spike back into lower triangle ====
//
            work(_(1,ns)) = V(1,_(1,ns));
            ElementType beta = work(1);
            ElementType tau;

            larfg(ns, beta, work(_(2,ns)), tau);
            work(1) = One;

            T(_(3,jw),_(1,jw-2)).lower() = Zero;

            const auto _v = work(_(1,ns));

            larf(Left,  _v, tau, T(_(1,ns),_(1,jw)), work(_(jw+1,jw+jw)));
            larf(Right, _v, tau, T(_(1,ns),_(1,ns)), work(_(jw+1,jw+ns)));
            larf(Right, _v, tau, V(_(1,jw),_(1,ns)), work(_(jw+1,jw+jw)));

            hrd(IndexType(1), ns, T_jw, work(_(1,jw-1)), work(_(jw+1,lWork)));
        }
//
//      ==== Copy updated reduced window into place ====
//
        if (kwTop>1) {
            H(kwTop,kwTop-1) = s*V(1,1);
        }
        H_kw.upper()  = T_jw.upper();
        H_kw.diag(-1) = T_jw.diag(-1);
//
//      ==== Accumulate orthogonal matrix in order update
//      .    H and Z, if requested.  ====
//
        if (ns>1 && s!=Zero) {
            ormhr(Right, NoTrans,
                  IndexType(1), ns,
                  T(_(1,ns),_(1,ns)),
                  work(_(1,ns-1)),
                  V(_(1,jw),_(1,ns)),
                  work(_(jw+1,lWork)));
        }
//
//      ==== Update vertical slab in H ====
//
        const IndexType lTop = (wantT) ? 1 : kTop;

        for (IndexType kRow=lTop; kRow<=kwTop-1; kRow+=nv) {
            const IndexType kLn = min(nv,kwTop-kRow);
            blas::mm(NoTrans, NoTrans, One,
                     H(_(kRow,kRow+kLn-1),_(kwTop,kBot)),
                     V(_(1,jw),_(1,jw)),
                     Zero,
                     WV(_(1,kLn),_(1,jw)));
            H(_(kRow,kRow+kLn-1),_(kwTop,kBot)) = WV(_(1,kLn),_(1,jw));
        }
//
//      ==== Update horizontal slab in H ====
//
        if (wantT) {
            for (IndexType kCol=kBot+1; kCol<=n; kCol+=nh) {
                const IndexType kLn = min(nh,n-kCol+1);
                blas::mm(ConjTrans, NoTrans, One,
                         V(_(1,jw),_(1,jw)),
                         H(_(kwTop,kBot),_(kCol,kCol+kLn-1)),
                         Zero,
                         T(_(1,jw),_(1,kLn)));
                H(_(kwTop,kBot),_(kCol,kCol+kLn-1)) = T(_(1,jw),_(1,kLn));
            }
        }
//
//      ==== Update vertical slab in Z ====
//
        if (wantZ) {
            for (IndexType kRow=iLoZ; kRow<=iHiZ; kRow+=nv) {
                const IndexType kLn = min(nv,iHiZ-kRow+1);
                blas::mm(NoTrans, NoTrans, One,
                         Z(_(kRow,kRow+kLn-1),_(kwTop,kBot)),
                         V(_(1,jw),_(1,jw)),
                         Zero,
                         WV(_(1,kLn),_(1,jw)));
                Z(_(kRow,kRow+kLn-1), _(kwTop,kBot)) = WV(_(1,kLn),_(1,jw));
            }
        }
    }
//
//  ==== Return the number of deflations ... ====
//
    nd = jw - ns;
//
//  ==== ... and the number of shifts. (Subtracting
//  .    INFQR from the spike length takes care
//  .    of the case of a rare QR failure while
//  .    calculating eigenvalues of the deflation
//  .    window.)  ====
//
    ns = ns - infoQR;
//
//   ==== Return optimal workspace. ====
//
    work(1) = lWorkOpt;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename MZ, typename VSH,
          typename MV, typename MT, typename MWV, typename VWORK>
void
laqr2_impl(bool                      wantT,
           bool                      wantZ,
           IndexType                 kTop,
           IndexType                 kBot,
           IndexType                 nw,
           GeMatrix<MH>              &H,
           IndexType                 iLoZ,
           IndexType                 iHiZ,
           GeMatrix<MZ>              &Z,
           IndexType                 &ns,
           IndexType                 &nd,
           DenseVector<VSH>          &sh,
           GeMatrix<MV>              &V,
           GeMatrix<MT>              &T,
           GeMatrix<MWV>             &WV,
           DenseVector<VWORK>        &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using cxxblas::abs1;

    typedef typename GeMatrix<MH>::ElementType                  ElementType;
    typedef typename ComplexTrait<ElementType>::PrimitiveType   PrimitiveType;

    const Underscore<IndexType> _;

    const ElementType    Zero(0), One(1);
    const PrimitiveType  RZero(0), ROne(1);
    const IndexType      n   = H.numRows();
    const IndexType      nh  = T.numCols();
    const IndexType      nv  = WV.numRows();
//
//  ==== Estimate optimal workspace. ====
//
    IndexType   jw       = min(nw, kBot-kTop+1);
    auto        T_jw     = T(_(1,jw),_(1,jw));
    auto        V_jw     = V(_(1,jw),_(1,jw));

    IndexType lWorkOpt;
    if (jw<=2) {
        lWorkOpt = 1;
    } else {
//
//      ==== Workspace query call to ZGEHRD ====
//
        IndexType lWork1 = hrd_wsq(IndexType(1), jw-1, T_jw);
//
//      ==== Workspace query call to ZUNMHR ====
//
        IndexType lWork2 = unmhr_wsq(Right, NoTrans, IndexType(1), jw-1, V_jw);
//
//      ==== Optimal workspace ====
//
        lWorkOpt = jw+max(lWork1,lWork2);
    }
//
//  ==== Apply worksize query ====
//
    if (work.length()==0) {
        work.resize(lWorkOpt);
    }
    const IndexType lWork = work.length();

//
//  ==== Nothing to do ...
//  ... for an empty active block ... ====
    ns = 0;
    nd = 0;
    work(1) = One;
    if (kTop>kBot) {
        return;
    }
//  ... nor for an empty deflation window. ====
    if (nw<1) {
        return;
    }
//
//  ==== Machine constants ====
//
    PrimitiveType safeMin   = lamch<PrimitiveType>(SafeMin);
    PrimitiveType safeMax   = ROne / safeMin;
    labad(safeMin, safeMax);
    const PrimitiveType ulp      = lamch<PrimitiveType>(Precision);
    const PrimitiveType smallNum = safeMin*(PrimitiveType(n)/ulp);
//
//  ==== Setup deflation window ====
//
    IndexType   kwTop   = kBot - jw + 1;
    auto        sh_kw   = sh(_(kwTop,kBot));
    auto        H_kw    = H(_(kwTop,kBot),_(kwTop,kBot));
    ElementType s;

    if( kwTop==kTop) {
        s = Zero;
    } else {
        s = H(kwTop, kwTop-1);
    }

    if (kBot==kwTop) {
//
//      ==== 1-by-1 deflation window: not much to do ====
//
        sh(kwTop) = H(kwTop,kwTop);
        ns = 1;
        nd = 0;
        if (abs1(s) <= max(smallNum, ulp*abs1(H(kwTop,kwTop)))) {
            ns = 0;
            nd = 1;
            if (kwTop > kTop) {
                H(kwTop,kwTop-1) = Zero;
            }
        }
        work(1) = One;
        return;
    }
//
//  ==== Convert to spike-triangular form.  (In case of a
//  .    rare QR failure, this routine continues to do
//  .    aggressive early deflation using that part of
//  .    the deflation window that converged using INFQR
//  .    here and there to keep track.) ====
//
    T_jw.upper()  = H_kw.upper();
    T_jw.diag(-1) = H_kw.diag(-1);

    V_jw          = Zero;
    V_jw.diag(0)  = One;

    IndexType infoQR;
    infoQR = lahqr(true, true,
                   IndexType(1), jw, T_jw,
                   sh_kw,
                   IndexType(1), jw, V_jw);
//
//  ==== Deflation detection loop ====
//
    ns = jw;
    IndexType iFirst;
    IndexType iLast = infoQR + 1;
    for (IndexType knt=infoQR+1; knt<=jw; ++knt) {
//
//      ==== Small spike tip deflation test ====
//
        PrimitiveType foo = abs1(T(ns,ns));
        if (foo==RZero) {
            foo = abs1(s);
        }
        if (abs1(s)*abs1(V(1,ns)) <= max(smallNum,ulp*foo)) {
//
//          ==== One more converged eigenvalue ====
//
            --ns;
        } else {
//
//          ==== One undeflatable eigenvalue.  Move it up out of the
//          .    way.   (ZTREXC can not fail in this case.) ====
//
            iFirst = ns;
            trexc(true, T_jw, V_jw, iFirst, iLast);
            ++iLast;
        }
    }
//
//      ==== Return to Hessenberg form ====
//
    if (ns==0) {
        s = Zero;
    }

    if (ns<jw) {
//
//      ==== sorting the diagonal of T improves accuracy for
//      .    graded matrices.  ====
//
        for (IndexType i=infoQR+1; i<=ns; ++i) {
            iFirst = i;
            for (IndexType j=i+1; j<=ns; ++j) {
                if (abs1(T(j,j))>abs1(T(iFirst,iFirst))) {
                    iFirst = j;
                }
            }
            iLast = i;
            if (iFirst != iLast) {
                trexc(true, T_jw, V_jw, iFirst, iLast);
            }
        }
    }
//
//  ==== Restore shift/eigenvalue array from T ====
//
    for (IndexType i=infoQR+1; i<=jw; ++i) {
        sh(kwTop+i-1) = T(i,i);
    }

    if (ns < jw || s==Zero) {

        if (ns > 1 && s != Zero) {
//
//          ==== Reflect spike back into lower triangle ====
//
            work(_(1,ns)) = V(1,_(1,ns));
            for (IndexType i=1; i<=ns; ++i) {
                work(i) = conj(work(i));
            }

            ElementType beta = work(1);
            ElementType tau;

            larfg(ns, beta, work(_(2,ns)), tau);
            work(1) = One;

            T(_(3,jw),_(1,jw-2)).lower() = Zero;

            const auto _v = work(_(1,ns));

            larf(Left,  _v, conj(tau), T(_(1,ns),_(1,jw)), work(_(jw+1,jw+jw)));
            larf(Right, _v, tau, T(_(1,ns),_(1,ns)), work(_(jw+1,jw+ns)));
            larf(Right, _v, tau, V(_(1,jw),_(1,ns)), work(_(jw+1,jw+jw)));

            hrd(IndexType(1), ns, T_jw, work(_(1,jw-1)), work(_(jw+1,lWork)));
        }
//
//      ==== Copy updated reduced window into place ====
//
        if (kwTop>1) {
            H(kwTop,kwTop-1) = s*conj(V(1,1));
        }
        H_kw.upper()  = T_jw.upper();
        H_kw.diag(-1) = T_jw.diag(-1);
//
//      ==== Accumulate orthogonal matrix in order update
//      .    H and Z, if requested.  ====
//
        if (ns>1 && s!=Zero) {
            unmhr(Right, NoTrans,
                  IndexType(1), ns,
                  T(_(1,ns),_(1,ns)),
                  work(_(1,ns-1)),
                  V(_(1,jw),_(1,ns)),
                  work(_(jw+1,lWork)));
        }
//
//      ==== Update vertical slab in H ====
//
        const IndexType lTop = (wantT) ? 1 : kTop;

        for (IndexType kRow=lTop; kRow<=kwTop-1; kRow+=nv) {
            const IndexType kLn = min(nv,kwTop-kRow);
            blas::mm(NoTrans, NoTrans, One,
                     H(_(kRow,kRow+kLn-1),_(kwTop,kBot)),
                     V(_(1,jw),_(1,jw)),
                     Zero,
                     WV(_(1,kLn),_(1,jw)));
            H(_(kRow,kRow+kLn-1),_(kwTop,kBot)) = WV(_(1,kLn),_(1,jw));
        }
//
//      ==== Update horizontal slab in H ====
//
        if (wantT) {
            for (IndexType kCol=kBot+1; kCol<=n; kCol+=nh) {
                const IndexType kLn = min(nh,n-kCol+1);
                blas::mm(ConjTrans, NoTrans, One,
                         V(_(1,jw),_(1,jw)),
                         H(_(kwTop,kBot),_(kCol,kCol+kLn-1)),
                         Zero,
                         T(_(1,jw),_(1,kLn)));
                H(_(kwTop,kBot),_(kCol,kCol+kLn-1)) = T(_(1,jw),_(1,kLn));
            }
        }
//
//      ==== Update vertical slab in Z ====
//
        if (wantZ) {
            for (IndexType kRow=iLoZ; kRow<=iHiZ; kRow+=nv) {
                const IndexType kLn = min(nv,iHiZ-kRow+1);
                blas::mm(NoTrans, NoTrans, One,
                         Z(_(kRow,kRow+kLn-1),_(kwTop,kBot)),
                         V(_(1,jw),_(1,jw)),
                         Zero,
                         WV(_(1,kLn),_(1,jw)));
                Z(_(kRow,kRow+kLn-1), _(kwTop,kBot)) = WV(_(1,kLn),_(1,jw));
            }
        }
    }
//
//  ==== Return the number of deflations ... ====
//
    nd = jw - ns;
//
//  ==== ... and the number of shifts. (Subtracting
//  .    INFQR from the spike length takes care
//  .    of the case of a rare QR failure while
//  .    calculating eigenvalues of the deflation
//  .    window.)  ====
//
    ns = ns - infoQR;

    //std::cerr << "work = " << work << std::endl;
    work(1) = lWorkOpt;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename IndexType, typename MH, typename MZ, typename VSR,
          typename VSI, typename MV, typename MT, typename MWV, typename VWORK>
void
laqr2_impl(bool                      wantT,
           bool                      wantZ,
           IndexType                 kTop,
           IndexType                 kBot,
           IndexType                 nw,
           GeMatrix<MH>              &H,
           IndexType                 iLoZ,
           IndexType                 iHiZ,
           GeMatrix<MZ>              &Z,
           IndexType                 &ns,
           IndexType                 &nd,
           DenseVector<VSR>          &sr,
           DenseVector<VSI>          &si,
           GeMatrix<MV>              &V,
           GeMatrix<MT>              &T,
           GeMatrix<MWV>             &WV,
           DenseVector<VWORK>        &work)
{
    cxxlapack::laqr2<IndexType>(wantT,
                                wantZ,
                                H.numRows(),
                                kTop,
                                kBot,
                                nw,
                                H.data(),
                                H.leadingDimension(),
                                iLoZ,
                                iHiZ,
                                Z.data(),
                                Z.leadingDimension(),
                                ns,
                                nd,
                                sr.data(),
                                si.data(),
                                V.data(),
                                V.leadingDimension(),
                                T.numCols(),
                                T.data(),
                                T.leadingDimension(),
                                WV.numRows(),
                                WV.data(),
                                WV.leadingDimension(),
                                work.data(),
                                work.length());
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename MZ, typename VSH,
          typename MV, typename MT, typename MWV, typename VWORK>
void
laqr2_impl(bool                      wantT,
           bool                      wantZ,
           IndexType                 kTop,
           IndexType                 kBot,
           IndexType                 nw,
           GeMatrix<MH>              &H,
           IndexType                 iLoZ,
           IndexType                 iHiZ,
           GeMatrix<MZ>              &Z,
           IndexType                 &ns,
           IndexType                 &nd,
           DenseVector<VSH>          &sh,
           GeMatrix<MV>              &V,
           GeMatrix<MT>              &T,
           GeMatrix<MWV>             &WV,
           DenseVector<VWORK>        &work)
{
    cxxlapack::laqr2<IndexType>(wantT,
                                wantZ,
                                H.numRows(),
                                kTop,
                                kBot,
                                nw,
                                H.data(),
                                H.leadingDimension(),
                                iLoZ,
                                iHiZ,
                                Z.data(),
                                Z.leadingDimension(),
                                ns,
                                nd,
                                sh.data(),
                                V.data(),
                                V.leadingDimension(),
                                T.numCols(),
                                T.data(),
                                T.leadingDimension(),
                                WV.numRows(),
                                WV.data(),
                                WV.leadingDimension(),
                                work.data(),
                                work.length());
}


//
//  Workspace query (real variant)
//
template <typename IndexType, typename MT>
typename RestrictTo<IsReal<typename MT::ElementType>::value,
         IndexType>::Type
laqr2_wsq_impl(IndexType                 kTop,
               IndexType                 kBot,
               IndexType                 nw,
               const GeMatrix<MT>        &T)
{
    typedef typename GeMatrix<MT>::ElementType  ElementType;

    ElementType     DUMMY, WORK;
    IndexType       NS, ND;
    IndexType       LWORK = -1;


    cxxlapack::laqr2<IndexType>(false,
                                false,
                                IndexType(1),
                                kTop,
                                kBot,
                                nw,
                                &DUMMY,                     // H.data()
                                IndexType(1),               // LDH
                                IndexType(0),               // iLoZ,
                                IndexType(0),               // iHiZ
                                &DUMMY,                     // H.data()
                                IndexType(1),               // LDZ
                                NS,
                                ND,
                                &DUMMY,                     // sr.data()
                                &DUMMY,                     // si.data()
                                &DUMMY,                     // V.data()
                                nw,                         // LDV
                                T.numCols(),
                                &DUMMY,
                                T.leadingDimension(),
                                nw,                         // NV
                                &DUMMY,                     // WV.data()
                                nw,                         // LDWV
                                &WORK,
                                LWORK);
    return WORK;
}

//
//  Workspace query (complex variant)
//
template <typename IndexType, typename MT>
typename RestrictTo<IsComplex<typename MT::ElementType>::value,
         IndexType>::Type
laqr2_wsq_impl(IndexType                 kTop,
               IndexType                 kBot,
               IndexType                 nw,
               const GeMatrix<MT>        &T)
{
    typedef typename GeMatrix<MT>::ElementType  ElementType;

    ElementType     DUMMY, WORK;
    IndexType       NS, ND;
    IndexType       LWORK = -1;


    cxxlapack::laqr2<IndexType>(false,
                                false,
                                IndexType(1),
                                kTop,
                                kBot,
                                nw,
                                &DUMMY,                     // H.data()
                                IndexType(1),               // LDH
                                IndexType(0),               // iLoZ,
                                IndexType(0),               // iHiZ
                                &DUMMY,                     // H.data()
                                IndexType(1),               // LDZ
                                NS,
                                ND,
                                &DUMMY,                     // sh.data()
                                &DUMMY,                     // V.data()
                                nw,                         // LDV
                                T.numCols(),
                                &DUMMY,
                                T.leadingDimension(),
                                nw,                         // NV
                                &DUMMY,                     // WV.data()
                                nw,                         // LDWV
                                &WORK,
                                LWORK);
    return cxxblas::real(WORK);
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
//
//  Real variant
//
template <typename IndexType, typename MH, typename MZ, typename VSR,
          typename VSI, typename MV, typename MT, typename MWV, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MH>::value
                 && IsRealGeMatrix<MZ>::value
                 && IsRealDenseVector<VSR>::value
                 && IsRealDenseVector<VSI>::value
                 && IsRealGeMatrix<MV>::value
                 && IsRealGeMatrix<MT>::value
                 && IsRealGeMatrix<MWV>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
laqr2(bool                      wantT,
      bool                      wantZ,
      IndexType                 kTop,
      IndexType                 kBot,
      IndexType                 nw,
      MH                        &&H,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      MZ                        &&Z,
      IndexType                 &ns,
      IndexType                 &nd,
      VSR                       &&sr,
      VSI                       &&si,
      MV                        &&V,
      MT                        &&T,
      MWV                       &&WV,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("laqr2");

    using std::max;

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type      MatrixH;
    typedef typename RemoveRef<MZ>::Type      MatrixZ;
    typedef typename RemoveRef<VSR>::Type     VectorSR;
    typedef typename RemoveRef<VSI>::Type     VectorSI;
    typedef typename RemoveRef<MV>::Type      MatrixV;
    typedef typename RemoveRef<MT>::Type      MatrixT;
    typedef typename RemoveRef<MWV>::Type     MatrixWV;
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    if (wantZ) {
        ASSERT(1<=iLoZ);
        ASSERT(iLoZ<=iHiZ);
        ASSERT(iHiZ<=n);

        ASSERT(Z.firstRow()==1);
        ASSERT(Z.firstCol()==1);
        ASSERT(Z.numRows()==n);
        ASSERT(Z.numCols()==n);
    }

    ASSERT(sr.length()==kBot);
    ASSERT(si.length()==kBot);

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);
    ASSERT(V.numRows()==nw);
    ASSERT(V.numCols()==nw);

    const IndexType nh = T.numCols();
    ASSERT(nh>=nw);

    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);

    const IndexType nv = WV.numRows();
    ASSERT(nv>=nw);

    ASSERT((work.length()==0) || (work.length()>=n));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixH::NoView        H_org      = H;
    typename MatrixZ::NoView        Z_org      = Z;
    IndexType                       ns_org     = ns;
    IndexType                       nd_org     = nd;
    typename VectorSR::NoView       sr_org     = sr;
    typename VectorSI::NoView       si_org     = si;
    typename MatrixV::NoView        V_org      = V;
    typename MatrixT::NoView        T_org      = T;
    typename MatrixWV::NoView       WV_org     = WV;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqr2_impl(wantT, wantZ, kTop, kBot, nw, H,
                              iLoZ, iHiZ, Z, ns, nd, sr, si,
                              V, T, WV, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixH::NoView        H_generic      = H;
    typename MatrixZ::NoView        Z_generic      = Z;
    IndexType                       ns_generic     = ns;
    IndexType                       nd_generic     = nd;
    typename VectorSR::NoView       sr_generic     = sr;
    typename VectorSI::NoView       si_generic     = si;
    typename MatrixV::NoView        V_generic      = V;
    typename MatrixT::NoView        T_generic      = T;
    typename MatrixWV::NoView       WV_generic     = WV;
    typename VectorWork::NoView     work_generic   = work;

//
//  restore output arguments
//
    H      = H_org;
    Z      = Z_org;
    ns     = ns_org;
    nd     = nd_org;
    sr     = sr_org;
    si     = si_org;
    V      = V_org;
    T      = T_org;
    WV     = WV_org;
    work   = work_org;

//
//  Compare results
//
    external::laqr2_impl(wantT, wantZ, kTop, kBot, nw, H,
                         iLoZ, iHiZ, Z, ns, nd, sr, si,
                         V, T, WV, work);

    bool failed = false;
    if (! isIdentical(H_generic, H, " H_generic", "H")) {
        std::cerr << "CXXLAPACK: H_generic = " << H_generic << std::endl;
        std::cerr << "F77LAPACK: H = " << H << std::endl;
        std::cerr << "Original:  H_org = " << H_org << std::endl;
        failed = true;
    }

    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
        failed = true;
    }

    if (! isIdentical(ns_generic, ns, "ns_generic", "ns")) {
        std::cerr << "CXXLAPACK: ns_generic = " << ns_generic << std::endl;
        std::cerr << "F77LAPACK: ns = " << ns << std::endl;
        failed = true;
    }

    if (! isIdentical(nd_generic, nd, "nd_generic", "nd")) {
        std::cerr << "CXXLAPACK: nd_generic = " << nd_generic << std::endl;
        std::cerr << "F77LAPACK: nd = " << nd << std::endl;
        failed = true;
    }

    if (! isIdentical(sr_generic, sr, "sr_generic", "_sr")) {
        std::cerr << "CXXLAPACK: sr_generic = " << sr_generic << std::endl;
        std::cerr << "F77LAPACK: sr = " << sr << std::endl;
        failed = true;
    }

    if (! isIdentical(si_generic, si, "si_generic", "si")) {
        std::cerr << "CXXLAPACK: si_generic = " << si_generic << std::endl;
        std::cerr << "F77LAPACK: si = " << si << std::endl;
        failed = true;
    }

    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }

    if (! isIdentical(T_generic, T, "T_generic", "T")) {
        std::cerr << "CXXLAPACK: T_generic = " << T_generic << std::endl;
        std::cerr << "F77LAPACK: T = " << T << std::endl;
        failed = true;
    }

    if (! isIdentical(WV_generic, WV, "WV_generic", "_WV")) {
        std::cerr << "CXXLAPACK: WV_generic = " << WV_generic << std::endl;
        std::cerr << "F77LAPACK: WV = " << WV << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "_work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr2.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: laqr2.tcc" << std::endl;
    }
#   endif
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename MZ, typename VSH,
          typename MV, typename MT, typename MWV, typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MH>::value
                 && IsComplexGeMatrix<MZ>::value
                 && IsComplexDenseVector<VSH>::value
                 && IsComplexGeMatrix<MV>::value
                 && IsComplexGeMatrix<MT>::value
                 && IsComplexGeMatrix<MWV>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
laqr2(bool                      wantT,
      bool                      wantZ,
      IndexType                 kTop,
      IndexType                 kBot,
      IndexType                 nw,
      MH                        &&H,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      MZ                        &&Z,
      IndexType                 &ns,
      IndexType                 &nd,
      VSH                       &&sh,
      MV                        &&V,
      MT                        &&T,
      MWV                       &&WV,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("laqr2 (complex)");

    using std::max;

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type      MatrixH;
    typedef typename RemoveRef<MZ>::Type      MatrixZ;
    typedef typename RemoveRef<VSH>::Type     VectorSH;
    typedef typename RemoveRef<MV>::Type      MatrixV;
    typedef typename RemoveRef<MT>::Type      MatrixT;
    typedef typename RemoveRef<MWV>::Type     MatrixWV;
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    if (wantZ) {
        ASSERT(1<=iLoZ);
        ASSERT(iLoZ<=iHiZ);
        ASSERT(iHiZ<=n);

        ASSERT(Z.firstRow()==1);
        ASSERT(Z.firstCol()==1);
        ASSERT(Z.numRows()==n);
        ASSERT(Z.numCols()==n);
    }

    ASSERT(sh.length()==kBot);

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);
    ASSERT(V.numRows()==nw);
    ASSERT(V.numCols()==nw);

    const IndexType nh = T.numCols();
    ASSERT(nh>=nw);

    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);

    const IndexType nv = WV.numRows();
    ASSERT(nv>=nw);

    ASSERT((work.length()==0) || (work.length()>=n));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixH::NoView        H_org      = H;
    typename MatrixZ::NoView        Z_org      = Z;
    IndexType                       ns_org     = ns;
    IndexType                       nd_org     = nd;
    typename VectorSH::NoView       sh_org     = sh;
    typename MatrixV::NoView        V_org      = V;
    typename MatrixT::NoView        T_org      = T;
    typename MatrixWV::NoView       WV_org     = WV;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::laqr2_impl(wantT, wantZ, kTop, kBot, nw, H, iLoZ, iHiZ, Z,
                              ns, nd, sh, V, T, WV, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixH::NoView        H_generic      = H;
    typename MatrixZ::NoView        Z_generic      = Z;
    IndexType                       ns_generic     = ns;
    IndexType                       nd_generic     = nd;
    typename VectorSH::NoView       sh_generic     = sh;
    typename MatrixV::NoView        V_generic      = V;
    typename MatrixT::NoView        T_generic      = T;
    typename MatrixWV::NoView       WV_generic     = WV;
    typename VectorWork::NoView     work_generic   = work;

//
//  restore output arguments
//
    H      = H_org;
    Z      = Z_org;
    ns     = ns_org;
    nd     = nd_org;
    sh     = sh_org;
    V      = V_org;
    T      = T_org;
    WV     = WV_org;
    work   = work_org;

//
//  Compare results
//
    external::laqr2_impl(wantT, wantZ, kTop, kBot, nw, H, iLoZ, iHiZ, Z,
                         ns, nd, sh, V, T, WV, work);

    bool failed = false;
    if (! isIdentical(H_generic, H, " H_generic", "H")) {
        std::cerr << "CXXLAPACK: H_generic = " << H_generic << std::endl;
        std::cerr << "F77LAPACK: H = " << H << std::endl;
        std::cerr << "Original:  H_org = " << H_org << std::endl;
        failed = true;
    }

    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
        failed = true;
    }

    if (! isIdentical(ns_generic, ns, "ns_generic", "ns")) {
        std::cerr << "CXXLAPACK: ns_generic = " << ns_generic << std::endl;
        std::cerr << "F77LAPACK: ns = " << ns << std::endl;
        failed = true;
    }

    if (! isIdentical(nd_generic, nd, "nd_generic", "nd")) {
        std::cerr << "CXXLAPACK: nd_generic = " << nd_generic << std::endl;
        std::cerr << "F77LAPACK: nd = " << nd << std::endl;
        failed = true;
    }

    if (! isIdentical(sh_generic, sh, "sh_generic", "_sh")) {
        std::cerr << "CXXLAPACK: sh_generic = " << sh_generic << std::endl;
        std::cerr << "F77LAPACK: sh = " << sh << std::endl;
        failed = true;
    }

    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }

    if (! isIdentical(T_generic, T, "T_generic", "T")) {
        std::cerr << "CXXLAPACK: T_generic = " << T_generic << std::endl;
        std::cerr << "F77LAPACK: T = " << T << std::endl;
        failed = true;
    }

    if (! isIdentical(WV_generic, WV, "WV_generic", "_WV")) {
        std::cerr << "CXXLAPACK: WV_generic = " << WV_generic << std::endl;
        std::cerr << "F77LAPACK: WV = " << WV << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "_work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr2.tcc" << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: laqr2.tcc" << std::endl;
    }
#   endif
}

//
//  Workspace query (real/complex variant)
//
template <typename IndexType, typename MT>
IndexType
laqr2_wsq(IndexType                 kTop,
          IndexType                 kBot,
          IndexType                 nw,
          const GeMatrix<MT>        &T)
{
    using std::max;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::laqr2_wsq_impl(kTop, kBot, nw, T);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::laqr2_wsq_impl(kTop, kBot, nw, T);

    if (info!=_info) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        ASSERT(0);
    }
#   endif
    return info;

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQR2_TCC
