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
      SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
      SUBROUTINE ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAQR0_TCC
#define FLENS_LAPACK_LA_LAQR0_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//
//  Workspace query (real/complex variant)
//
template <typename IndexType, typename MH>
IndexType
laqr0_wsq_impl(bool                  wantT,
               bool                  wantZ,
               IndexType             iLo,
               IndexType             iHi,
               const GeMatrix<MH>    &H)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MH>::ElementType  T;

    const IndexType nTiny   = 11;
    const IndexType n       = H.numRows();

    if ((n==0) || (n<=nTiny)) {
        return 1;
    }

    char job[3];
    job[0] = (wantT) ? 'S' : 'E';
    job[1] = (wantZ) ? 'V' : 'N';
    job[2] = 0;
//
//  ==== NWR = recommended deflation window size.  At this
//  .    point,  N .GT. NTINY = 11, so there is enough
//  .    subdiagonal workspace for NWR.GE.2 as required.
//  .    (In fact, there is enough subdiagonal space for
//  .    NWR.GE.3.) ====
//
    IndexType nwr = ilaenv<T>(13, "LAQR0", job, n, iLo, iHi, -1);
    nwr = max(IndexType(2), nwr);
    nwr = min(min(IndexType(iHi-iLo+1), (n-1)/3), nwr);
//
//  ==== NSR = recommended number of simultaneous shifts.
//  .    At this point N .GT. NTINY = 11, so there is at
//  .    enough subdiagonal workspace for NSR to be even
//  .    and greater than or equal to two as required. ====
//
    IndexType nsr = ilaenv<T>(15, "LAQR0", job, n, iLo, iHi, -1);
    nsr = min(min(nsr, (n+6)/9), IndexType(iHi-iLo));
    nsr = max(IndexType(2), nsr-(nsr%2));
//
//  ==== Estimate optimal workspace ====
//
//  ==== Workspace query call to DLAQR3 ====
//
    IndexType lWorkOpt = laqr3_wsq(IndexType(iLo), IndexType(iHi), nwr+1, H);
//
//  ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
//
    return max(3*nsr/2, lWorkOpt);
}

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
laqr0_impl(bool                  wantT,
           bool                  wantZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VWR>      &wr,
           DenseVector<VWI>      &wi,
           IndexType             iLoZ,
           IndexType             iHiZ,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    using flens::min;
    using std::abs;
    using std::max;
    using std::swap;

    typedef typename GeMatrix<MH>::ElementType  T;

    const Underscore<IndexType>     _;

    const IndexType n       = H.numRows();

//  ==== Matrices of order NTINY or smaller must be processed by
//  .    DLAHQR because of insufficient subdiagonal scratch space.
//  .    (This is a hard limit.) ====
    const IndexType nTiny   = 11;

//  ==== Exceptional deflation windows:  try to cure rare
//  .    slow convergence by varying the size of the
//  .    deflation window after KEXNW iterations. ====
    const IndexType kexNw   = 5;

//
//  ==== Exceptional shifts: try to cure rare slow convergence
//  .    with ad-hoc exceptional shifts every KEXSH iterations.
//  .    ====
    const IndexType kexSh   = 6;

//
//  ==== The constants WILK1 and WILK2 are used to form the
//  .    exceptional shifts. ====
    const T         wilk1 = T(0.75),
                    wilk2 = T(-0.4375);

    const T         Zero(0), One(1);

    IndexType info          = 0;
    IndexType nDec          = -1;

    IndexType lWork;
    IndexType lWorkOpt;

//
//  ==== Perform and apply a workspace query if necessary ====
//
    if (work.length()==0) {
        lWorkOpt = laqr0_wsq(wantT, wantZ, iLo, iHi, H);
        work.resize(lWorkOpt);
    }

    lWork = work.length();

//
//  ==== Quick return for N = 0: nothing to do. ====
//
    if (n==0) {
        work(1) = One;
        return info;
    }

    if (n<=nTiny) {
//
//      ==== Tiny matrices must use DLAHQR. ====
//
        lWorkOpt = 1;
        info = lahqr(wantT, wantZ, iLo, iHi, H, wr, wi, iLoZ, iHiZ, Z);
    } else {
//
//      ==== Use small bulge multi-shift QR with aggressive early
//      .    deflation on larger-than-tiny matrices. ====
//
//      ==== Hope for the best. ====
//
        info = 0;
//
//      ==== Set up job flags for ILAENV. ====
//
        char job[3];
        job[0] = (wantT) ? 'S' : 'E';
        job[1] = (wantZ) ? 'V' : 'N';
        job[2] = 0;
//
//      ==== NWR = recommended deflation window size.  At this
//      .    point,  N .GT. NTINY = 11, so there is enough
//      .    subdiagonal workspace for NWR.GE.2 as required.
//      .    (In fact, there is enough subdiagonal space for
//      .    NWR.GE.3.) ====
//
        IndexType nwr = ilaenv<T>(13, "LAQR0", job, n, iLo, iHi, lWork);
        nwr = max(IndexType(2), nwr);
        nwr = min(min(IndexType(iHi-iLo+1), (n-1)/3), nwr);
//
//      ==== NSR = recommended number of simultaneous shifts.
//      .    At this point N .GT. NTINY = 11, so there is at
//      .    enough subdiagonal workspace for NSR to be even
//      .    and greater than or equal to two as required. ====
//
        IndexType nsr = ilaenv<T>(15, "LAQR0", job, n, iLo, iHi, lWork);
        nsr = min(min(nsr, (n+6)/9), IndexType(iHi-iLo));
        nsr = max(IndexType(2), nsr - (nsr%2));
//
//      ==== Estimate optimal workspace ====
//
//      ==== Workspace query call to DLAQR3 ====
//
        lWorkOpt = laqr3_wsq(iLo, iHi, nwr+1, H);
//
//      ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
//
        lWorkOpt = max(3*nsr/2, lWorkOpt);
//
//      ==== DLAHQR/DLAQR0 crossover point ====
//
        IndexType nMin = ilaenv<T>(12, "LAQR0", job, n, iLo, iHi, lWork);
        nMin = max(nTiny, nMin);
//
//      ==== Nibble crossover point ====
//
        IndexType nibble = ilaenv<T>(14, "LAQR0", job, n, iLo, iHi, lWork);
        nibble = max(IndexType(0), nibble);
//
//      ==== Accumulate reflections during ttswp?  Use block
//      .    2-by-2 structure during matrix-matrix multiply? ====
//
        IndexType kacc22 = ilaenv<T>(16, "LAQR0", job, n, iLo, iHi, lWork);
        kacc22 = max(IndexType(0), kacc22);
        kacc22 = min(IndexType(2), kacc22);
//
//      ==== NWMAX = the largest possible deflation window for
//      .    which there is sufficient workspace. ====
//
        IndexType nwMax = min((n-1)/3, lWork/2);
        IndexType nw = nwMax;
//
//      ==== NSMAX = the Largest number of simultaneous shifts
//      .    for which there is sufficient workspace. ====
//
        IndexType nsMax = min((n+6 )/9, 2*lWork/3);
        nsMax -= nsMax % 2;
//
//      ==== NDFL: an iteration count restarted at deflation. ====
//
        IndexType nDfl = 1;
//
//      ==== ITMAX = iteration limit ====
//
        IndexType itMax = max(IndexType(30), 2*kexSh)
                        * max(IndexType(10), iHi-iLo+1);
//
//      ==== Last row and column in the active block ====
//
        IndexType kBot = iHi;
//
//      ==== Main Loop ====
//
        IndexType it;
        for (it=1; it<=itMax; ++it) {
//
//          ==== Done when KBOT falls below ILO ====
//
            if (kBot<iLo) {
                break;
            }
//
//          ==== Locate active block ====
//
            IndexType k;
            for (k=kBot; k>=iLo+1; --k) {
                if (H(k,k-1)==Zero) {
                    break;
                }
            }
            ASSERT(k==iLo || H(k,k-1)==Zero);
            const IndexType kTop = k;
//
//          ==== Select deflation window size:
//          .    Typical Case:
//          .      If possible and advisable, nibble the entire
//          .      active block.  If not, use size MIN(NWR,NWMAX)
//          .      or MIN(NWR+1,NWMAX) depending upon which has
//          .      the smaller corresponding subdiagonal entry
//          .      (a heuristic).
//          .
//          .    Exceptional Case:
//          .      If there have been no deflations in KEXNW or
//          .      more iterations, then vary the deflation window
//          .      size.   At first, because, larger windows are,
//          .      in general, more powerful than smaller ones,
//          .      rapidly increase the window to the maximum possible.
//          .      Then, gradually reduce the window size. ====
//
            IndexType nh = kBot - kTop + 1;
            IndexType nwUpBd = min(nh, nwMax);
            if (nDfl<kexNw) {
                nw = min(nwUpBd, nwr);
            } else {
                nw = min(nwUpBd, 2*nw);
            }
            if (nw<nwMax) {
                if (nw>=nh-1) {
                    nw = nh;
                } else {
                    const IndexType kwTop = kBot - nw + 1;
                    if (abs(H(kwTop,kwTop-1))>abs(H(kwTop-1, kwTop-2))) {
                        ++nw;
                    }
                }
            }
            if (nDfl<kexNw) {
                nDec = -1;
            } else if (nDec>=0 || nw>=nwUpBd) {
                ++nDec;
                if (nw-nDec<2) {
                    nDec = 0;
                }
                nw -= nDec;
            }
//
//          ==== Aggressive early deflation:
//          .    split workspace under the subdiagonal into
//          .      - an nw-by-nw work array V in the lower
//          .        left-hand-corner,
//          .      - an NW-by-at-least-NW-but-more-is-better
//          .        (NW-by-NHO) horizontal work array along
//          .        the bottom edge,
//          .      - an at-least-NW-but-more-is-better (NHV-by-NW)
//          .        vertical work array along the left-hand-edge.
//          .        ====
//
            auto _V     = H(_(n-nw+1,    n), _(   1,     nw));
            auto _T     = H(_(n-nw+1,    n), _(nw+1, n-nw-1));
            auto _WV    = H(_(  nw+2, n-nw), _(   1,     nw));
//
//          ==== Aggressive early deflation ====
//
            IndexType   ls, ld;

            laqr3(wantT, wantZ, kTop, kBot, nw, H, iLoZ, iHiZ, Z, ls, ld,
                  wr(_(1,kBot)), wi(_(1,kBot)), _V, _T, _WV, work);
//
//          ==== Adjust KBOT accounting for new deflations. ====
//
            kBot -= ld;
//
//          ==== KS points to the shifts. ====
//
            IndexType ks = kBot - ls + 1;
//
//          ==== Skip an expensive QR sweep if there is a (partly
//          .    heuristic) reason to expect that many eigenvalues
//          .    will deflate without it.  Here, the QR sweep is
//          .    skipped if many eigenvalues have just been deflated
//          .    or if the remaining active block is small.
//
            if ((ld==0)
             || ((100*ld<=nw*nibble) && (kBot-kTop+1>min(nMin, nwMax))))
            {
//
//              ==== NS = nominal number of simultaneous shifts.
//              .    This may be lowered (slightly) if DLAQR3
//              .    did not provide that many shifts. ====
//
                IndexType ns = min(nsMax, nsr, max(IndexType(2), kBot-kTop));
                ns -= ns % 2;
//
//              ==== If there have been no deflations
//              .    in a multiple of KEXSH iterations,
//              .    then try exceptional shifts.
//              .    Otherwise use shifts provided by
//              .    DLAQR3 above or from the eigenvalues
//              .    of a trailing principal submatrix. ====
//
                if (nDfl%kexSh==0) {
                    ks = kBot - ns + 1;
                    for (IndexType i=kBot; i>=max(ks+1,kTop+2); i-=2) {
                        const T ss = abs(H(i,i-1)) + abs(H(i-1,i-2));
                        T aa = wilk1*ss + H(i,i);
                        T bb = ss;
                        T cc = wilk2*ss;
                        T dd = aa;
                        T cs, sn;
                        lanv2(aa, bb, cc, dd,
                              wr(i-1), wi(i-1), wr(i), wi(i),
                              cs, sn);
                    }
                    if (ks==kTop) {
                        wr(ks+1)    = H(ks+1, ks+1);
                        wi(ks+1)    = Zero;
                        wr(ks)      = wr(ks+1);
                        wi(ks)      = wi(ks+1);
                    }
                } else {
//
//                  ==== Got NS/2 or fewer shifts? Use DLAQR4 or
//                  .    DLAHQR on a trailing principal submatrix to
//                  .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
//                  .    there is enough space below the subdiagonal
//                  .    to fit an NS-by-NS scratch array.) ====
//
                    if (kBot-ks+1<=ns/2) {
                        ks = kBot - ns +1;
                        H(_(ks,kBot),_(1,ns)) = H(_(ks,kBot),_(ks,kBot));
                        if (ns>nMin) {
                            // TODO: avoid the need for ZDummy
                            typename GeMatrix<MZ>::NoView ZDummy;
                            ks += laqr4(false, false,
                                        IndexType(1), ns,
                                        H(_(ks,kBot),_(1,ns)),
                                        wr(_(ks,kBot)), wi(_(ks,kBot)),
                                        IndexType(1), IndexType(1),
                                        ZDummy, work);
                        } else {
                            // TODO: avoid the need for ZDummy
                            typename GeMatrix<MZ>::NoView ZDummy;
                            ks += lahqr(false, false,
                                        IndexType(1), ns,
                                        H(_(ks,kBot),_(1,ns)),
                                        wr(_(ks,kBot)), wi(_(ks,kBot)),
                                        IndexType(1), IndexType(1),
                                        ZDummy);
                        }
//
//                      ==== In case of a rare QR failure use
//                      .    eigenvalues of the trailing 2-by-2
//                      .    principal submatrix.  ====
//
                        if (ks>=kBot) {
                            T aa = H(kBot-1,kBot-1);
                            T cc = H(kBot,  kBot-1);
                            T bb = H(kBot-1,kBot);
                            T dd = H(kBot,  kBot);
                            T cs, sn;
                            lanv2(aa, bb, cc, dd,
                                  wr(kBot-1), wi(kBot-1), wr(kBot), wi(kBot),
                                  cs, sn);
                            ks = kBot - 1;
                        }
                    }

                    if (kBot-ks+1>ns) {
//
//                      ==== Sort the shifts (Helps a little)
//                      .    Bubble sort keeps complex conjugate
//                      .    pairs together. ====
//
                        bool sorted = false;
                        for (IndexType k=kBot; k>=ks+1; --k) {
                            if (sorted) {
                                break;
                            }
                            sorted = true;
                            for (IndexType i=ks; i<=k-1; ++i) {
                                if (abs(wr(i))+abs(wi(i))
                                    < abs(wr(i+1))+abs(wi(i+1)))
                                {
                                    sorted = false;
                                    swap(wr(i), wr(i+1));
                                    swap(wi(i), wi(i+1));
                                }
                            }
                        }
                    }
//
//                  ==== Shuffle shifts into pairs of real shifts
//                  .    and pairs of complex conjugate shifts
//                  .    assuming complex conjugate shifts are
//                  .    already adjacent to one another. (Yes,
//                  .    they are.)  ====
//
                    for (IndexType i=kBot; i>=ks+2; i-=2) {
                        if (wi(i)!=-wi(i-1)) {
                            T tmp  = wr(i);
                            wr(i)   = wr(i-1);
                            wr(i-1) = wr(i-2);
                            wr(i-2) = tmp;

                            tmp     = wi(i);
                            wi(i)   = wi(i-1);
                            wi(i-1) = wi(i-2);
                            wi(i-2) = tmp;
                        }
                    }
                }
//
//              ==== If there are only two shifts and both are
//              .    real, then use only one.  ====
//
                if (kBot-ks+1==2) {
                    if (wi(kBot)==0) {
                        const T _H = H(kBot,kBot);
                        if (abs(wr(kBot)-_H) < abs(wr(kBot-1)-_H)) {
                            wr(kBot-1) = wr(kBot);
                        } else {
                            wr(kBot) = wr(kBot-1);
                        }
                    }
                }
//
//              ==== Use up to NS of the the smallest magnatiude
//              .    shifts.  If there aren't NS shifts available,
//              .    then use them all, possibly dropping one to
//              .    make the number of shifts even. ====
//
                ns = min(ns, kBot-ks+1);
                ns -= ns % 2;
                ks = kBot - ns + 1;
//
//              ==== Small-bulge multi-shift QR sweep:
//              .    split workspace under the subdiagonal into
//              .    - a KDU-by-KDU work array U in the lower
//              .      left-hand-corner,
//              .    - a KDU-by-at-least-KDU-but-more-is-better
//              .      (KDU-by-NHo) horizontal work array WH along
//              .      the bottom edge,
//              .    - and an at-least-KDU-but-more-is-better-by-KDU
//              .      (NVE-by-KDU) vertical work WV arrow along
//              .      the left-hand-edge. ====
//
                IndexType kdu   = 3*ns - 3;
                IndexType ku    = n - kdu + 1;
                IndexType kwv   = kdu + 4;
                IndexType nho   = (n-kdu+1-4) - (kdu+1) + 1;

                typedef typename GeMatrix<MH>::View GeMatrixView;
                GeMatrixView _V(IndexType(3), ns/2, work(_(1,3*ns/2)));
                auto _U     = H(_( ku,    n), _(    1,    kdu));
                auto _WV    = H(_(kwv,n-kdu), _(    1,    kdu));
                auto _WH    = H(_( ku,    n), _(kdu+1,kdu+nho));
//
//              ==== Small-bulge multi-shift QR sweep ====
//
                laqr5(wantT, wantZ, kacc22, kTop, kBot, ns,
                      wr(_(ks,kBot)), wi(_(ks,kBot)), H,
                      iLoZ, iHiZ, Z, _V, _U, _WV, _WH);
            }
//
//          ==== Note progress (or the lack of it). ====
//
            if (ld>0) {
                nDfl = 1;
            } else {
                ++nDfl;
            }
//
//          ==== End of main loop ====
//
        }
//
//      ==== Iteration limit exceeded.  Set INFO to show where
//      .    the problem occurred and exit. ====
//
        if (it>itMax) {
            info = kBot;
        }
    }

    work(1) = lWorkOpt;
    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
IndexType
laqr0_impl(bool                  wantT,
           bool                  wantZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VW>       &w,
           IndexType             iLoZ,
           IndexType             iHiZ,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    using std::abs;
    using std::max;
    using std::swap;
    using cxxblas::abs1;
    using flens::min;

    typedef typename GeMatrix<MH>::ElementType       T;
    typedef typename ComplexTrait<T>::PrimitiveType  PT;

    const Underscore<IndexType>     _;

    const IndexType n       = H.numRows();
//
//  ==== Matrices of order NTINY or smaller must be processed by
//  .    ZLAHQR because of insufficient subdiagonal scratch space.
//  .    (This is a hard limit.) ====
    const IndexType     nTiny = 11;
//
//  ==== Exceptional deflation windows:  try to cure rare
//  .    slow convergence by varying the size of the
//  .    deflation window after KEXNW iterations. ====
    const IndexType     kexNw = 5;
//
//     ==== Exceptional shifts: try to cure rare slow convergence
//     .    with ad-hoc exceptional shifts every KEXSH iterations.
//     .    ====
    const IndexType     kexSh = 6;
//
//     ==== The constant WILK1 is used to form the exceptional
//     .    shifts. ====
    const PT            wilk1 = PT(0.75);

    const T             Zero(0), One(1);
    const PT            Two(2);

    IndexType info          = 0;

    IndexType lWork;
    IndexType lWorkOpt;

    IndexType nDec          = -1;
//
//  ==== Perform and apply a workspace query if necessary ====
//
    if (work.length()==0) {
        lWorkOpt = laqr0_wsq(wantT, wantZ, iLo, iHi, H);
        work.resize(lWorkOpt);
    }

    lWork = work.length();

//
//  ==== Quick return for N = 0: nothing to do. ====
//
    if (n==0) {
        work(1) = One;
        return info;
    }

    if (n<=nTiny) {
//
//      ==== Tiny matrices must use ZLAHQR. ====
//
        lWorkOpt = 1;
        info = lahqr(wantT, wantZ, iLo, iHi, H, w, iLoZ, iHiZ, Z);
    } else {
//
//      ==== Use small bulge multi-shift QR with aggressive early
//      .    deflation on larger-than-tiny matrices. ====
//
//      ==== Hope for the best. ====
//
        info = 0;
//
//      ==== Set up job flags for ILAENV. ====
//
        char job[3];
        job[0] = (wantT) ? 'S' : 'E';
        job[1] = (wantZ) ? 'V' : 'N';
        job[2] = 0;
//
//      ==== NWR = recommended deflation window size.  At this
//      .    point,  N .GT. NTINY = 11, so there is enough
//      .    subdiagonal workspace for NWR.GE.2 as required.
//      .    (In fact, there is enough subdiagonal space for
//      .    NWR.GE.3.) ====
//
        IndexType nwr = ilaenv<T>(13, "LAQR0", job, n, iLo, iHi, lWork);
        nwr = max(IndexType(2), nwr);
        nwr = min(iHi-iLo+1, (n-1)/3, nwr);
//
//       ==== NSR = recommended number of simultaneous shifts.
//       .    At this point N .GT. NTINY = 11, so there is at
//       .    enough subdiagonal workspace for NSR to be even
//       .    and greater than or equal to two as required. ====
//
        IndexType nsr = ilaenv<T>(15, "LAQR0", job, n, iLo, iHi, lWork);
        nsr = min(nsr, (n+6)/9, iHi-iLo);
        nsr = max(IndexType(2), nsr - (nsr%2));
//
//      ==== Estimate optimal workspace ====
//
//      ==== Workspace query call to ZLAQR3 ====
//
        lWorkOpt = laqr3_wsq(iLo, iHi, nwr+1, H);
//
//      ==== Optimal workspace = MAX(ZLAQR5, ZLAQR3) ====
//
        lWorkOpt = max(3*nsr/2, lWorkOpt);
//
//      ==== ZLAHQR/ZLAQR0 crossover point ====
//
        IndexType nMin = ilaenv<T>(12, "laqr0", job, n, iLo, iHi, lWork);
        nMin = max(nTiny, nMin);
//
//      ==== Nibble crossover point ====
//
        IndexType nibble = ilaenv<T>(14, "laqr0", job, n, iLo, iHi, lWork);
        nibble = max(IndexType(0), nibble);
//
//      ==== Accumulate reflections during ttswp?  Use block
//      .    2-by-2 structure during matrix-matrix multiply? ====
//
        IndexType kacc22 = ilaenv<T>(16, "laqr0", job, n, iLo, iHi, lWork);
        kacc22 = max(IndexType(0), kacc22);
        kacc22 = min(IndexType(2), kacc22);
//
//      ==== NWMAX = the largest possible deflation window for
//      .    which there is sufficient workspace. ====
//
        IndexType nwMax = min((n-1)/3, lWork/2);
        IndexType nw = nwMax;
//
//      ==== NSMAX = the Largest number of simultaneous shifts
//      .    for which there is sufficient workspace. ====
//
        IndexType nsMax = min((n+6)/9, 2*lWork/3);
        nsMax -= nsMax % 2;
//
//      ==== NDFL: an iteration count restarted at deflation. ====
//
        IndexType nDfl = 1;
//
//      ==== ITMAX = iteration limit ====
//
        IndexType itMax = max(IndexType(30), 2*kexSh)
                        * max(IndexType(10), iHi-iLo+1);
//
//      ==== Last row and column in the active block ====
//
        IndexType kBot = iHi;
//
//      ==== Main Loop ====
//
        IndexType it;
        for (it=1; it<=itMax; ++it) {
//
//          ==== Done when KBOT falls below ILO ====
//
            if (kBot<iLo) {
                break;
            }
//
//          ==== Locate active block ====
//
            IndexType k;
            for (k=kBot; k>=iLo+1; --k) {
                if (H(k,k-1)==Zero) {
                    break;
                }
            }
            ASSERT(k==iLo || H(k,k-1)==Zero);
            const IndexType kTop = k;
//
//          ==== Select deflation window size:
//          .    Typical Case:
//          .      If possible and advisable, nibble the entire
//          .      active block.  If not, use size MIN(NWR,NWMAX)
//          .      or MIN(NWR+1,NWMAX) depending upon which has
//          .      the smaller corresponding subdiagonal entry
//          .      (a heuristic).
//          .
//          .    Exceptional Case:
//          .      If there have been no deflations in KEXNW or
//          .      more iterations, then vary the deflation window
//          .      size.   At first, because, larger windows are,
//          .      in general, more powerful than smaller ones,
//          .      rapidly increase the window to the maximum possible.
//          .      Then, gradually reduce the window size. ====
//
            IndexType nh = kBot - kTop + 1;
            IndexType nwUpBd = min(nh, nwMax);
            if (nDfl<kexNw) {
                nw = min(nwUpBd, nwr);
            } else {
                nw = min(nwUpBd, 2*nw);
            }
            if (nw<nwMax) {
                if (nw>=nh-1) {
                    nw = nh;
                } else {
                    const IndexType kwTop = kBot - nw + 1;
                    if (abs1(H(kwTop,kwTop-1)) > abs1(H(kwTop-1,kwTop-2))) {
                        ++nw;
                    }
                }
            }
            if (nDfl<kexNw) {
                nDec = -1;
            } else if (nDec>=0 || nw>=nwUpBd) {
                ++nDec;
                if (nw-nDec<2) {
                    nDec = 0;
                }
                nw -= nDec;
            }
//
//          ==== Aggressive early deflation:
//          .    split workspace under the subdiagonal into
//          .      - an nw-by-nw work array V in the lower
//          .        left-hand-corner,
//          .      - an NW-by-at-least-NW-but-more-is-better
//          .        (NW-by-NHO) horizontal work array along
//          .        the bottom edge,
//          .      - an at-least-NW-but-more-is-better (NHV-by-NW)
//          .        vertical work array along the left-hand-edge.
//          .        ====
//
            auto _V     = H(_(n-nw+1,    n), _(   1,     nw));
            auto _T     = H(_(n-nw+1,    n), _(nw+1, n-nw-1));
            auto _WV    = H(_(  nw+2, n-nw), _(   1,     nw));
//
//          ==== Aggressive early deflation ====
//
            IndexType   ls, ld;

            laqr3(wantT, wantZ, kTop, kBot, nw, H, iLoZ, iHiZ, Z, ls, ld,
                  w(_(1,kBot)), _V, _T, _WV, work);
//
//          ==== Adjust KBOT accounting for new deflations. ====
//
            kBot -= ld;
//
//          ==== KS points to the shifts. ====
//
            IndexType ks = kBot - ls + 1;
//
//          ==== Skip an expensive QR sweep if there is a (partly
//          .    heuristic) reason to expect that many eigenvalues
//          .    will deflate without it.  Here, the QR sweep is
//          .    skipped if many eigenvalues have just been deflated
//          .    or if the remaining active block is small.
//
            if ((ld==0)
             || ((100*ld <= nw*nibble) && (kBot-kTop+1 > min(nMin, nwMax))))
            {
//
//              ==== NS = nominal number of simultaneous shifts.
//              .    This may be lowered (slightly) if ZLAQR3
//              .    did not provide that many shifts. ====
//
                IndexType ns = min(nsMax, nsr, max(IndexType(2), kBot-kTop));
                ns -= ns % 2;
//
//              ==== If there have been no deflations
//              .    in a multiple of KEXSH iterations,
//              .    then try exceptional shifts.
//              .    Otherwise use shifts provided by
//              .    ZLAQR3 above or from the eigenvalues
//              .    of a trailing principal submatrix. ====
//
                if (nDfl%kexSh==0) {
                    ks = kBot - ns + 1;
                    for (IndexType i=kBot; i>=ks+1; i-=2) {
                        w(i)   = H(i,i) + wilk1*abs1(H(i,i-1));
                        w(i-1) = w(i);
                    }
                } else {
//
//                  ==== Got NS/2 or fewer shifts? Use ZLAQR4 or
//                  .    ZLAHQR on a trailing principal submatrix to
//                  .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
//                  .    there is enough space below the subdiagonal
//                  .    to fit an NS-by-NS scratch array.) ====
//
                    if (kBot-ks+1<=ns/2) {
                        ks = kBot - ns +1;
                        H(_(ks,kBot),_(1,ns)) = H(_(ks,kBot),_(ks,kBot));

                        if (ns>nMin) {
                            // TODO: avoid the need for ZDummy
                            typename GeMatrix<MZ>::NoView ZDummy;
                            ks += laqr4(false, false,
                                        IndexType(1), ns,
                                        H(_(ks,kBot),_(1,ns)),
                                        w(_(ks,kBot)),
                                        IndexType(1), IndexType(1),
                                        ZDummy, work);
                        } else {
                            // TODO: avoid the need for ZDummy
                            typename GeMatrix<MZ>::NoView ZDummy;
                            ks += lahqr(false, false,
                                        IndexType(1), ns,
                                        H(_(ks,kBot),_(1,ns)),
                                        w(_(ks,kBot)),
                                        IndexType(1), IndexType(1),
                                        ZDummy);
                        }
//                      ==== In case of a rare QR failure use
//                      .    eigenvalues of the trailing 2-by-2
//                      .    principal submatrix.  Scale to avoid
//                      .    overflows, underflows and subnormals.
//                      .    (The scale factor S can not be zero,
//                      .    because H(KBOT,KBOT-1) is nonzero.) ====
//
                        if (ks>=kBot) {
                            PT s = abs1(H(kBot-1, kBot-1)) +
                                   abs1(H(kBot,   kBot-1)) +
                                   abs1(H(kBot-1, kBot  )) +
                                   abs1(H(kBot,   kBot  ));
                            T aa = H(kBot-1, kBot-1) / s;
                            T cc = H(kBot,   kBot-1) / s;
                            T bb = H(kBot-1, kBot  ) / s;
                            T dd = H(kBot,   kBot  ) / s;

                            T tr2 = (aa+dd) / Two;
                            T det = (aa-tr2) * (dd-tr2) - bb*cc;
                            T rtDisc = sqrt(-det);

                            w(kBot-1) = (tr2+rtDisc)*s;
                            w(kBot)   = (tr2-rtDisc)*s;

                            ks = kBot - 1;
                        }
                    }

                    if (kBot-ks+1>ns) {
//
//                      ==== Sort the shifts (Helps a little) ====
//
                        bool sorted = false;
                        for (IndexType k=kBot; k>=ks+1; --k) {
                            if (sorted) {
                                break;
                            }
                            sorted = true;
                            for (IndexType i=ks; i<=k-1; ++i) {
                                if (abs1(w(i))<abs1(w(i+1))) {
                                    sorted = false;
                                    swap(w(i), w(i+1));
                                }
                            }
                        }
                    }
                }
//
//              ==== If there are only two shifts, then use
//              .    only one.  ====
//
                if (kBot-ks+1==2) {
                    if (abs1(w(kBot)-H(kBot,kBot))
                      < abs1(w(kBot-1)-H(kBot,kBot)))
                    {
                        w(kBot-1) = w(kBot);
                    } else {
                        w(kBot) = w(kBot-1);
                    }
                }
//
//              ==== Use up to NS of the the smallest magnatiude
//              .    shifts.  If there aren't NS shifts available,
//              .    then use them all, possibly dropping one to
//              .    make the number of shifts even. ====
//
                ns = min(ns, kBot-ks+1);
                ns -= ns % 2;
                ks = kBot - ns + 1;

//
//              ==== Small-bulge multi-shift QR sweep:
//              .    split workspace under the subdiagonal into
//              .    - a KDU-by-KDU work array U in the lower
//              .      left-hand-corner,
//              .    - a KDU-by-at-least-KDU-but-more-is-better
//              .      (KDU-by-NHo) horizontal work array WH along
//              .      the bottom edge,
//              .    - and an at-least-KDU-but-more-is-better-by-KDU
//              .      (NVE-by-KDU) vertical work WV arrow along
//              .      the left-hand-edge. ====
//
                IndexType kdu   = 3*ns - 3;
                IndexType ku    = n - kdu + 1;
                IndexType kwv   = kdu + 4;
                IndexType nho   = (n-kdu+1-4) - (kdu+1) + 1;

                typedef typename GeMatrix<MH>::View GeMatrixView;
                GeMatrixView _V(IndexType(3), ns/2, work(_(1,3*ns/2)));
                auto _U     = H(_( ku,    n), _(    1,    kdu));
                auto _WV    = H(_(kwv,n-kdu), _(    1,    kdu));
                auto _WH    = H(_( ku,    n), _(kdu+1,kdu+nho));
//
//              ==== Small-bulge multi-shift QR sweep ====
//
                laqr5(wantT, wantZ, kacc22, kTop, kBot, ns,
                      w(_(ks,kBot)), H, iLoZ, iHiZ, Z,
                      _V, _U, _WV, _WH);

            }
//
//          ==== Note progress (or the lack of it). ====
//
            if (ld>0) {
                nDfl = 1;
            } else {
                ++nDfl;
            }
//
//          ==== End of main loop ====
//
        }
//
//      ==== Iteration limit exceeded.  Set INFO to show where
//      .    the problem occurred and exit. ====
//
        if (it>itMax) {
            info = kBot;
        }
    }

    work(1) = T(lWorkOpt, 0);
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//
//  Real variant
//
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
laqr0_impl(bool                  wantT,
           bool                  wantZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VWR>      &wr,
           DenseVector<VWI>      &wi,
           IndexType             iLoZ,
           IndexType             iHiZ,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    IndexType  info;
    info = cxxlapack::laqr0<IndexType>(wantT,
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
                                       Z.leadingDimension(),
                                       work.data(),
                                       work.length());
    ASSERT(info>=0);
    return info;
}

//
//  Workspace query (real variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsReal<typename MH::ElementType>::value,
         IndexType>::Type
laqr0_wsq_impl(bool                  wantT,
               bool                  wantZ,
               IndexType             iLo,
               IndexType             iHi,
               const GeMatrix<MH>    &H)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    T                   WORK, DUMMY;
    const IndexType     LWORK   = -1;

    cxxlapack::laqr0<IndexType>(wantT,
                                wantZ,
                                H.numRows(),
                                iLo,
                                iHi,
                                &DUMMY,
                                H.leadingDimension(),
                                &DUMMY,
                                &DUMMY,
                                IndexType(1),
                                IndexType(1),
                                &DUMMY,
                                IndexType(1),
                                &WORK,
                                LWORK);
    return WORK;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
IndexType
laqr0_impl(bool                  wantT,
           bool                  wantZ,
           IndexType             iLo,
           IndexType             iHi,
           GeMatrix<MH>          &H,
           DenseVector<VW>       &w,
           IndexType             iLoZ,
           IndexType             iHiZ,
           GeMatrix<MZ>          &Z,
           DenseVector<VWORK>    &work)
{
    IndexType  info;
    info = cxxlapack::laqr0<IndexType>(wantT,
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
                                       Z.leadingDimension(),
                                       work.data(),
                                       work.length());
    ASSERT(info>=0);
    return info;
}

//
//  Workspace query (complex variant)
//
template <typename IndexType, typename MH>
typename RestrictTo<IsComplex<typename MH::ElementType>::value,
         IndexType>::Type
laqr0_wsq_impl(bool                  wantT,
               bool                  wantZ,
               IndexType             iLo,
               IndexType             iHi,
               const GeMatrix<MH>    &H)
{
    typedef typename GeMatrix<MH>::ElementType  T;

    T                   WORK, DUMMY;
    const IndexType     LWORK   = -1;

    cxxlapack::laqr0<IndexType>(wantT,
                                wantZ,
                                H.numRows(),
                                iLo,
                                iHi,
                                &DUMMY,
                                H.leadingDimension(),
                                &DUMMY,
                                IndexType(1),
                                IndexType(1),
                                &DUMMY,
                                IndexType(1),
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
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
typename RestrictTo<IsRealGeMatrix<MH>::value
                 && IsRealDenseVector<VWR>::value
                 && IsRealDenseVector<VWI>::value
                 && IsRealGeMatrix<MZ>::value
                 && IsRealDenseVector<VWORK>::value,
         IndexType>::Type
laqr0(bool                      wantT,
      bool                      wantZ,
      IndexType                 iLo,
      IndexType                 iHi,
      MH                        &&H,
      VWR                       &&wr,
      VWI                       &&wi,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      MZ                        &&Z,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("laqr0");

    using std::max;
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type        MatrixH;
    typedef typename RemoveRef<VWR>::Type       VectorWR;
    typedef typename RemoveRef<VWI>::Type       VectorWI;
    typedef typename RemoveRef<MZ>::Type        MatrixZ;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    if (n>0) {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    } else {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    }

    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()>=iHi);

    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()>=iHi);

    ASSERT(1<=iLoZ);
    ASSERT(iLoZ<=iLo);
    ASSERT(iHi<=iHiZ);
    ASSERT(iHiZ<=n);

    ASSERT(Z.firstRow()==1);
    ASSERT(Z.firstCol()==1);
    ASSERT(Z.numRows()>=iHi);
    ASSERT(Z.numCols()>=iHi);

    ASSERT((work.length()==0) || (work.length()>=n));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixH::NoView       H_org       = H;
    typename VectorWR::NoView      wr_org      = wr;
    typename VectorWI::NoView      wi_org      = wi;
    typename MatrixZ::NoView       Z_org       = Z;
    typename VectorWork::NoView    work_org    = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::laqr0_impl(wantT, wantZ, iLo, iHi, H,
                                               wr, wi, iLoZ, iHiZ, Z, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixH::NoView       H_generic       = H;
    typename VectorWR::NoView      wr_generic      = wr;
    typename VectorWI::NoView      wi_generic      = wi;
    typename MatrixZ::NoView       Z_generic       = Z;
    typename VectorWork::NoView    work_generic    = work;
//
//  restore output arguments
//
    H    = H_org;
    wr   = wr_org;
    wi   = wi_org;
    Z    = Z_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::laqr0_impl(wantT, wantZ, iLo, iHi, H, wr, wi,
                                           iLoZ, iHiZ, Z, work);

    bool failed = false;
    if (! isIdentical(H_generic, H, "H_generic", "H")) {
        std::cerr << "CXXLAPACK: H_generic = " << H_generic << std::endl;
        std::cerr << "F77LAPACK: H = " << H << std::endl;
        failed = true;
    }

    if (! isIdentical(wr_generic, wr, "wr_generic", "wr")) {
        std::cerr << "CXXLAPACK: wr_generic = " << wr_generic << std::endl;
        std::cerr << "F77LAPACK: wr = " << wr << std::endl;
        failed = true;
    }

    if (! isIdentical(wi_generic, wi, "wi_generic", "wi")) {
        std::cerr << "CXXLAPACK: wi_generic = " << wi_generic << std::endl;
        std::cerr << "F77LAPACK: wi = " << wi << std::endl;
        failed = true;
    }

    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr0.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: laqr0.tcc" << std::endl;
    }
#   endif

    return info;
}

//
//  Complex variant
//
template <typename IndexType, typename MH, typename VW, typename MZ,
          typename VWORK>
typename RestrictTo<IsComplexGeMatrix<MH>::value
                 && IsComplexDenseVector<VW>::value
                 && IsComplexGeMatrix<MZ>::value
                 && IsComplexDenseVector<VWORK>::value,
         IndexType>::Type
laqr0(bool                      wantT,
      bool                      wantZ,
      IndexType                 iLo,
      IndexType                 iHi,
      MH                        &&H,
      VW                        &&w,
      IndexType                 iLoZ,
      IndexType                 iHiZ,
      MZ                        &&Z,
      VWORK                     &&work)
{
    LAPACK_DEBUG_OUT("laqr0");

    using std::max;
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MH>::Type        MatrixH;
    typedef typename RemoveRef<VW>::Type        VectorW;
    typedef typename RemoveRef<MZ>::Type        MatrixZ;
    typedef typename RemoveRef<VWORK>::Type     VectorWork;
#   endif
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    if (n>0) {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    } else {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    }

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()>=iHi);

    ASSERT(1<=iLoZ);
    ASSERT(iLoZ<=iLo);
    ASSERT(iHi<=iHiZ);
    ASSERT(iHiZ<=n);

    ASSERT(Z.firstRow()==1);
    ASSERT(Z.firstCol()==1);
    ASSERT(Z.numRows()>=iHi);
    ASSERT(Z.numCols()>=iHi);

    ASSERT((work.length()==0) || (work.length()>=n));
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixH::NoView       H_org       = H;
    typename VectorW::NoView       w_org       = w;
    typename MatrixZ::NoView       Z_org       = Z;
    typename VectorWork::NoView    work_org    = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::laqr0_impl(wantT, wantZ, iLo, iHi, H,
                                               w, iLoZ, iHiZ, Z, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixH::NoView       H_generic       = H;
    typename VectorW::NoView       w_generic       = w;
    typename MatrixZ::NoView       Z_generic       = Z;
    typename VectorWork::NoView    work_generic    = work;
//
//  restore output arguments
//
    H    = H_org;
    w    = w_org;
    Z    = Z_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::laqr0_impl(wantT, wantZ, iLo, iHi, H, w,
                                           iLoZ, iHiZ, Z, work);

    bool failed = false;
    if (! isIdentical(H_generic, H, "H_generic", "H")) {
        std::cerr << "CXXLAPACK: H_generic = " << H_generic << std::endl;
        std::cerr << "F77LAPACK: H = " << H << std::endl;
        failed = true;
    }

    if (! isIdentical(w_generic, w, "w_generic", "w")) {
        std::cerr << "CXXLAPACK: w_generic = " << w_generic << std::endl;
        std::cerr << "F77LAPACK: w = " << w << std::endl;
        failed = true;
    }

    if (! isIdentical(Z_generic, Z, "Z_generic", "Z")) {
        std::cerr << "CXXLAPACK: Z_generic = " << Z_generic << std::endl;
        std::cerr << "F77LAPACK: Z = " << Z << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: laqr0.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: laqr0.tcc" << std::endl;
    }
#   endif

    return info;
}


//
//  Workspace query (real/complex variant)
//
template <typename IndexType, typename MH>
IndexType
laqr0_wsq(bool                  wantT,
          bool                  wantZ,
          IndexType             iLo,
          IndexType             iHi,
          const GeMatrix<MH>    &H)
{
    using std::max;
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());

    const IndexType n = H.numRows();

    if (n>0) {
        ASSERT(1<=iLo);
        ASSERT(iLo<=iHi);
        ASSERT(iHi<=n);
    } else {
        ASSERT(iLo==1);
        ASSERT(iHi==0);
    }
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::laqr0_wsq_impl(wantT, wantZ, iLo, iHi, H);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::laqr0_wsq_impl(wantT, wantZ, iLo, iHi, H);

    if (info!=_info) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        ASSERT(0);
    }
#   endif
    return info;
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAQR0_TCC
