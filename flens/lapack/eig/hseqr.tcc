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
      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
     $                   LDZ, WORK, LWORK, INFO )
 *
 *  -- LAPACK computational routine (version 3.2.2) --
 *     Univ. of Tennessee, Univ. of California Berkeley,
 *     Univ. of Colorado Denver and NAG Ltd..
 *     June 2010
 */

#ifndef FLENS_LAPACK_EIG_HSEQR_TCC
#define FLENS_LAPACK_EIG_HSEQR_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

using namespace HSEQR;

using std::max;
using std::min;

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
hseqr(Job job, ComputeZ computeZ, IndexType iLo, IndexType iHi, MH &&H,
      VWR &&wr, VWI &&wi, MZ &&Z, VWORK &&work)
{
    hseqr(job, computeZ, iLo, iHi, H, wr, wi, Z, work);
}

//-- hseqr ---------------------------------------------------------------------
template <typename IndexType, typename MH, typename VWR, typename VWI,
          typename MZ, typename VWORK>
IndexType
hseqr(Job job, ComputeZ computeZ, IndexType iLo, IndexType iHi, GeMatrix<MH> &H,
      DenseVector<VWR> &wr, DenseVector<VWI> &wi, GeMatrix<MZ> &Z,
      DenseVector<VWORK> &work)
{
    ASSERT(H.firstRow()==1);
    ASSERT(H.firstCol()==1);
    ASSERT(H.numRows()==H.numCols());
    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==H.numCols());
    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==H.numCols());
    ASSERT((computeZ==No) || (Z.numRows()==H.numCols()));
    ASSERT((computeZ==No) || (Z.numCols()==H.numCols()));
    ASSERT((work.length()==0) || (work.length()>=H.numCols()));

    typedef typename GeMatrix<MH>::ElementType   T;
    typedef typename GeMatrix<MH>::View          HLMatrix;
    typedef typename GeMatrix<MH>::VectorView    WorkL;

    const IndexType nTiny = 11;
//  ==== NL allocates some local workspace to help small matrices
//  .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
//  .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
//  .    mended.  (The default value of NMIN is 75.)  Using NL = 49
//  .    allows up to six simultaneous shifts and a 16-by-16
//  .    deflation window.  ====
    const       IndexType nl = 49;
    T           hlBuffer[nl*nl], worklBuffer[nl];
    HLMatrix    Hl = typename HLMatrix::Engine(nl, nl, hlBuffer, nl);
    WorkL       workl = typename WorkL::Engine(nl, worklBuffer);

    const Underscore<IndexType>     _;
    const IndexType                 n = H.numCols();
    const T                         Zero(0), One(1);

    IndexType info = 0;

//
//  ==== Decode and check the input parameters. ====
//
    const bool wantT = (job==Schur);
    const bool initZ = (computeZ==Init);
    const bool wantZ = (computeZ!=No);

    if (work.length()>0) {
        work(1) = n;
    } else {
//
//      ==== Apply a workspace query ====
//
        // TODO: uncomment this
        //laqr0(wantT, wantZ, iLo, iHi, H, wr, wi, iLo, iHi, Z, work);
        ASSERT(0);
    }

    if (n==0) {
//
//      ==== Quick return in case n = 0; nothing to do. ====
//
        return info;
    }
//
//  ==== copy eigenvalues isolated by DGEBAL ====
//
    for (IndexType i=1; i<iLo; ++i) {
        wr(i) = H(i,i);
        wi(i) = Zero;
    }

    for (IndexType i=iHi+1; i<=n; ++i) {
        wr(i) = H(i,i);
        wi(i) = Zero;
    }
//
//  ==== Initialize Z, if requested ====
//
    if (initZ) {
        // TODO: implement laset and uncomment this
        // laset(Zero, One, Z);
        Z = Zero;
        Z.diag(0) = One;
    }
//
//  ==== Quick return if possible ====
//
    if (iLo==iHi) {
        wr(iLo) = H(iLo, iLo);
        wi(iLo) = Zero;
        return info;
    }
//
//  ==== DLAHQR/DLAQR0 crossover point ====
//
    char opt[3];
    opt[0] = (job==Eigenvalues) ? 'E' : 'S';
    if (computeZ==No) {
        opt[1] = 'N';
    } else if (computeZ==Init) {
        opt[1] = 'I';
    } else {
        opt[1] = 'V';
    }
    opt[2] = 0;
    // TODO: uncomment this
    //IndexType nMin = ilaenv<T>(12, "HSEQR", opt, n, iLo, iHi, work.length());
    IndexType nMin = 20;
    
    nMin= max(nTiny, nMin);
//
//  ==== DLAQR0 for big matrices; DLAHQR for small ones ====
//
    if (n>nMin) {
        // TODO: uncomment this
        //info = laqr0(wantT, wantZ, iLo, iHi, H, wr, wi, iLo, iHi, Z, work);
        ASSERT(0);
    } else {
//
//      ==== Small matrix ====
//
        info = lahqr(wantT, wantZ, iLo, iHi, H, wr, wi, iLo, iHi, Z);

        if (info>0) {
//
//          ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
//          .    when DLAHQR fails. ====
//
            const IndexType kBot = info;
            if (n>=nl) {
//
//              ==== Larger matrices have enough subdiagonal scratch
//              .    space to call DLAQR0 directly. ====
//
                // TODO: uncomment this
                /*
                info = laqr0(wantT, wantZ,
                             iLo, kBot, H, wr, wi,
                             iLo, iHi, Z, work)
                */
                ASSERT(0);
            } else {
//
//              ==== Tiny matrices don't have enough subdiagonal
//              .    scratch space to benefit from DLAQR0.  Hence,
//              .    tiny matrices must be copied into a larger
//              .    array before calling DLAQR0. ====
//
                blas::copy(NoTrans, H, Hl(_(1,n),_(1,n)));
                Hl(n+1, n) = Zero;
                Hl(_(1,nl),_(n+1,nl)) = Zero;
                // TODO: uncomment this
                /*
                info = laqr0(wantT, wantZ,
                             iLo, kBot, Hl, wr, wi,
                             iLo, iHi, Z, workl);
                */
                if (wantT || (info==0)) {
                    blas::copy(NoTrans, Hl(_(1,n),_(1,n)), H);
                }
            }
        }
    }
//
//  ==== Clear out the trash, if necessary. ====
//
    if ((wantT || (info==0)) && (n>2)) {
        H(_(3,n-1),_(1,n-1)).lower() = Zero;
    }
//
//  ==== Ensure reported workspace size is backward-compatible with
//  .    previous LAPACK versions. ====
//
    work(1) = max(T(max(IndexType(1),n)), work(1));
    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_HSEQR_TCC
