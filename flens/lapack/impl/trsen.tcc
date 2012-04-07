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
       SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,
      $                   M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_IMPL_TRSEN_TCC
#define FLENS_LAPACK_IMPL_TRSEN_TCC 1

#include <cmath>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
namespace generic {

template <typename SELECT, typename MT, typename IndexType>
Pair<IndexType>
trsen_wsq_impl(TRSEN::Job                job,
               const DenseVector<SELECT> &select,
               GeMatrix<MT>              &T,
               IndexType                 &m)
{
    using std::max;

    typedef typename GeMatrix<MT>::ElementType  ElementType;
    const ElementType   Zero(0);

    const IndexType n = T.numRows();

//
//  Decode and test the input parameters
//
    const bool wantBH = (job==TRSEN::Both);
    const bool wantSP = (job==TRSEN::InvariantSubspaceOnly) || wantBH;
//
//  Set M to the dimension of the specified invariant subspace,
//  and test LWORK and LIWORK.
//
    m = 0;
    bool pair = false;
    for (IndexType k=1; k<=n; ++k) {
        if (pair) {
            pair = false;
        } else {
            if (k<n) {
                if (T(k+1,k)==Zero) {
                    if (select(k)) {
                        ++m;
                    }
                } else {
                    pair = true;
                    if (select(k) || select(k+1)) {
                        m += 2;
                    }
                }
            } else {
                if (select(n)) {
                    ++m;
                }
            }
        }
    }

    IndexType n1 = m;
    IndexType n2 = n - m;
    IndexType nn = n1*n2;

    IndexType workMin = 0,
              iWorkMin = 0;

    if (wantSP) {
        workMin  = max(IndexType(1), 2*nn);
        iWorkMin = max(IndexType(1), nn);
    } else if (job==TRSEN::None) {
        workMin  = max(IndexType(1), n);
        iWorkMin = 1;
    } else if (job==TRSEN::EigenvaluesOnly) {
        workMin  = max(IndexType(1), nn);
        iWorkMin = 1;
    }

    return Pair<IndexType>(workMin, iWorkMin);
}

template <typename SELECT, typename MT, typename MQ, typename WR, typename WI,
          typename IndexType, typename S, typename SEP,
          typename WORK, typename IWORK>
IndexType
trsen_impl(TRSEN::Job                job,
           bool                      computeQ,
           const DenseVector<SELECT> &select,
           GeMatrix<MT>              &T,
           GeMatrix<MQ>              &Q,
           DenseVector<WR>           &wr,
           DenseVector<WI>           &wi,
           IndexType                 &m,
           S                         &s,
           SEP                       &sep,
           DenseVector<WORK>         &work,
           DenseVector<IWORK>        &iwork)
{
    using std::abs;

    typedef typename GeMatrix<MT>::ElementType  ElementType;
    const ElementType   Zero(0), One(1);

    const IndexType n = T.numRows();
    const Underscore<IndexType> _;

    IndexType info = 0;

//
//  .. Local Arrays ..
//  this array is used to save variables between calls to lacn2
//
    IndexType _isaveData[3] = {0, 0, 0};
    DenseVectorView<IndexType>
        isave = typename DenseVectorView<IndexType>::Engine(3, _isaveData);
//
//  Decode and test the input parameters
//
    const bool wantBH = (job==TRSEN::Both);
    const bool wantS  = (job==TRSEN::EigenvaluesOnly) || wantBH;
    const bool wantSP = (job==TRSEN::InvariantSubspaceOnly) || wantBH;
//
//  Set M to the dimension of the specified invariant subspace,
//  and compute minimal size of WORK and IWORK
//
    auto wsq = trsen_wsq_impl(job, select, T, m);
    const IndexType workMin  = wsq.first;
    const IndexType iWorkMin = wsq.second;

    IndexType lWork  = work.length();
    IndexType liWork = iwork.length();

    if (lWork<workMin && lWork>0) {
        ASSERT(0);
        return -1;
    } else if (liWork<iWorkMin && liWork>0) {
        ASSERT(0);
        return -2;
    }

    if (lWork==0) {
        work.resize(workMin);
    }
    if (liWork==0) {
        iwork.resize(iWorkMin);
    }

    work(1)  = workMin;
    iwork(1) = iWorkMin;

    IndexType n1 = m;
    IndexType n2 = n - m;
    IndexType nn = n1*n2;

    bool pair;
    IndexType ks;
    ElementType scale;

//
//  Quick return if possible.
//
    if (m==n || m==0) {
        if (wantS) {
            s = One;
        }
        if (wantSP) {
            sep = lan(OneNorm, T);
        }
        goto DONE;
    }
//
//  Collect the selected blocks at the top-left corner of T.
//
    ks = 0;
    pair = false;
    for (IndexType k=1; k<=n; ++k) {
        if (pair) {
            pair = false;
        } else {
            bool swap = select(k);
            if (k<n) {
                if (T(k+1,k)!=Zero) {
                    pair = true;
                    swap = swap || select(k+1);
                }
            }
            if (swap) {
                ++ks;
//
//              Swap the K-th block to position KS.
//
                IndexType iErr = 0;
                IndexType kk = k;
                if (k!=ks) {
                    trexc(computeQ, T, Q, kk, ks, work(_(1,n)));
                }
                if (iErr==1 || iErr==2) {
//
//                  Blocks too close to swap: exit.
//
                    info = 1;
                    if (wantS) {
                        s = Zero;
                    }
                    if (wantSP) {
                        sep = Zero;
                    }
                    goto DONE;
                }
                if (pair) {
                    ++ks;
                }
            }
        }
    }

    if (wantS) {
        GeMatrixView<ElementType>    Work(n1, n2, work, n1);
//
//      Solve Sylvester equation for R:
//
//      T11*R - R*T22 = scale*T12
//
        Work = T(_(1,n1),_(n1+1,n1+n2));
        auto T11 = T(_(1,n1),_(1,n1));
        auto T22 = T(_(n1+1,n1+n2),_(n1+1,n1+n2));

        trsyl(NoTrans, NoTrans, -1, T11, T22, Work, scale);
//
//      Estimate the reciprocal of the condition number of the cluster
//      of eigenvalues.
//
        ElementType rNorm = lan(FrobeniusNorm, Work);
        if (rNorm==Zero) {
            s = One;
        } else {
            s = scale / (sqrt(scale*scale/rNorm + rNorm)*sqrt(rNorm));
        }
    }

    if (wantSP) {
//
//      Estimate sep(T11,T22).
//
        ElementType est = Zero;
        IndexType   kase = 0;

        do {
            auto _v     = work(_(nn+1,nn+nn));
            auto _x     = work(_(1,nn));
            auto _isgn  = iwork(_(1,nn));

            lacn2(_v, _x, _isgn, est, kase, isave);
            if (kase==0) {
                break;
            }

            GeMatrixView<ElementType>   Work(n1, n2, work, n1);
            auto T11 = T(_(1,n1),_(1,n1));
            auto T22 = T(_(n1+1,n1+n2),_(n1+1,n1+n2));

            if (kase==1) {
//
//              Solve  T11*R - R*T22 = scale*X.
//
                trsyl(NoTrans, NoTrans, -1, T11, T22, Work, scale);
            } else {
//
//              Solve T11**T*R - R*T22**T = scale*X.
//
                trsyl(Trans, Trans, -1, T11, T22, Work, scale);
            }
        } while (true);

        sep = scale / est;
    }

    DONE:

//
//  Store the output eigenvalues in WR and WI.
//
    for (IndexType k=1; k<=n; ++k) {
        wr(k) = T(k,k);
        wi(k) = Zero;
    }
    for (IndexType k=1; k<=n-1; ++k) {
        if (T(k+1,k)!=Zero) {
            wi(k)   = sqrt(abs(T(k,k+1)))*sqrt(abs(T(k+1,k)));
            wi(k+1) = -wi(k);
        }
    }

    work(1) = workMin;
    iwork(1) = iWorkMin;
    return info;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename SELECT, typename MT, typename MQ, typename WR, typename WI,
          typename IndexType, typename S, typename SEP,
          typename WORK, typename IWORK>
IndexType
trsen_impl(TRSEN::Job                job,
           bool                      computeQ,
           const DenseVector<SELECT> &select,
           GeMatrix<MT>              &T,
           GeMatrix<MQ>              &Q,
           DenseVector<WR>           &wr,
           DenseVector<WI>           &wi,
           IndexType                 &m,
           S                         &s,
           SEP                       &sep,
           DenseVector<WORK>         &work,
           DenseVector<IWORK>        &iwork)
{
    IndexType  info;
    info = cxxlapack::trsen<IndexType>(getF77Char(job),
                                       computeQ ? 'V' : 'N',
                                       select.data(),
                                       T.numRows(),
                                       T.data(),
                                       T.leadingDimension(),
                                       Q.data(),
                                       Q.leadingDimension(),
                                       wr.data(),
                                       wi.data(),
                                       m,
                                       s,
                                       sep,
                                       work.data(),
                                       work.length(),
                                       iwork.data(),
                                       iwork.length());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename SELECT, typename MT, typename MQ, typename WR, typename WI,
          typename IndexType, typename S, typename SEP,
          typename WORK, typename IWORK>
IndexType
trsen(TRSEN::Job                job,
      bool                      computeQ,
      const DenseVector<SELECT> &select,
      GeMatrix<MT>              &T,
      GeMatrix<MQ>              &Q,
      DenseVector<WR>           &wr,
      DenseVector<WI>           &wi,
      IndexType                 &m,
      S                         &s,
      SEP                       &sep,
      DenseVector<WORK>         &work,
      DenseVector<IWORK>        &iwork)
{
    LAPACK_DEBUG_OUT("trsen");
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());

    const IndexType n = T.numRows();
    if (computeQ) {
        ASSERT(Q.firstRow()==1);
        ASSERT(Q.firstCol()==1);
        ASSERT(Q.numRows()==n);
        ASSERT(Q.numCols()==n);
    }

    ASSERT(wr.firstIndex()==1);
    ASSERT(wr.length()==n);

    ASSERT(wi.firstIndex()==1);
    ASSERT(wi.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    const typename GeMatrix<MT>::NoView           T_org     = T;
    const typename GeMatrix<MQ>::NoView           Q_org     = Q;
    const typename DenseVector<WR>::NoView        wr_org    = wr;
    const typename DenseVector<WI>::NoView        wi_org    = wi;
    const IndexType                               m_org     = m;
    const S                                       s_org     = s;
    const SEP                                     sep_org   = sep;
    const typename DenseVector<WORK>::NoView      work_org  = work;
    const typename DenseVector<IWORK>::NoView     iwork_org = iwork;
#   endif

//
//  Call implementation
//
    IndexType info =  LAPACK_SELECT::trsen_impl(job, computeQ, select,
                                                T, Q, wr, wi,
                                                m, s, sep, work, iwork);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by generic implementation
//
    const typename GeMatrix<MT>::NoView           T_generic     = T;
    const typename GeMatrix<MQ>::NoView           Q_generic     = Q;
    const typename DenseVector<WR>::NoView        wr_generic    = wr;
    const typename DenseVector<WI>::NoView        wi_generic    = wi;
    const IndexType                               m_generic     = m;
    const S                                       s_generic     = s;
    const SEP                                     sep_generic   = sep;
    const typename DenseVector<WORK>::NoView      work_generic  = work;
    const typename DenseVector<IWORK>::NoView     iwork_generic = iwork;

//
//  restore output arguments
//
    T     = T_org;
    Q     = Q_org;
    wr    = wr_org;
    wi    = wi_org;
    m     = m_org;
    s     = s_org;
    sep   = sep_org;
    work  = work_org;
    iwork = iwork_org;

//
//  Compare generic results with results from the native implementation
//
    IndexType _info =  external::trsen_impl(job, computeQ, select,
                                            T, Q, wr, wi,
                                            m, s, sep, work, iwork);
    bool failed = false;
    if (! isIdentical(T_generic, T, "T_generic", "T")) {
        std::cerr << "CXXLAPACK: T_generic = " << T_generic << std::endl;
        std::cerr << "F77LAPACK: T = " << T << std::endl;
        failed = true;
    }
    if (! isIdentical(Q_generic, Q, "Q_generic", "Q")) {
        std::cerr << "CXXLAPACK: Q_generic = " << Q_generic << std::endl;
        std::cerr << "F77LAPACK: Q = " << Q << std::endl;
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
    if (! isIdentical(m_generic, m, "m_generic", "m")) {
        std::cerr << "CXXLAPACK: m_generic = " << m_generic << std::endl;
        std::cerr << "F77LAPACK: m = " << m << std::endl;
        failed = true;
    }
    if (! isIdentical(s_generic, s, "s_generic", "s")) {
        std::cerr << "CXXLAPACK: s_generic = " << s_generic << std::endl;
        std::cerr << "F77LAPACK: s = " << s << std::endl;
        failed = true;
    }
    if (! isIdentical(sep_generic, sep, "sep_generic", "sep")) {
        std::cerr << "CXXLAPACK: sep_generic = " << sep_generic << std::endl;
        std::cerr << "F77LAPACK: sep = " << sep << std::endl;
        failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }
    if (! isIdentical(iwork_generic, iwork, "iwork_generic", "iwork")) {
        std::cerr << "CXXLAPACK: iwork_generic = "
                  << iwork_generic << std::endl;
        std::cerr << "F77LAPACK: iwork = " << iwork << std::endl;
        failed = true;
    }
    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info =  " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: trsen.tcc" << std::endl;
        ASSERT(0);
    } else {
//      std::cerr << "passed: trsen.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------

template <typename SELECT, typename MT, typename MQ, typename WR, typename WI,
          typename IndexType, typename S, typename SEP,
          typename WORK, typename IWORK>
IndexType
trsen(TRSEN::Job                job,
      bool                      computeQ,
      const DenseVector<SELECT> &select,
      MT                        &&T,
      MQ                        &&Q,
      WR                        &&wr,
      WI                        &&wi,
      IndexType                 &&m,
      S                         &&s,
      SEP                       &&sep,
      WORK                      &&work,
      IWORK                     &&iwork)
{
    trsen(job, computeQ, select, T, Q, wr, wi, m, s, sep, work, iwork);
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRSEN_TCC
