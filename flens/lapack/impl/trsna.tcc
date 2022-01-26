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
      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,
     $                   INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_TRSNA_TCC
#define FLENS_LAPACK_IMPL_TRSNA_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
namespace generic {

template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename M, typename MM,
          typename MWORK, typename VIWORK>
void
trsna_impl(TRSNA::Job                    job,
           TRSNA::HowMany                howMany,
           const DenseVector<VSELECT>    &select,
           const GeMatrix<MT>            &T,
           const GeMatrix<MVL>           &VL_,
           const GeMatrix<MVR>           &VR_,
           DenseVector<VS>               &s,
           DenseVector<VSEP>             &sep,
           const MM                      &mm,
           M                             &m,
           GeMatrix<MWORK>               &Work,
           DenseVector<VIWORK>           &iWork)
{
    using std::abs;

    typedef typename GeMatrix<MT>::ElementType  ElementType;
    typedef typename GeMatrix<MT>::IndexType    IndexType;

    const Underscore<IndexType> _;

    const IndexType n = T.numRows();

    const ElementType Zero(0), One(1), Two(2);
//
//  Local Arrays
//
    IndexType iSaveData[3] = {0, 0, 0};
    DenseVectorView<IndexType>
        iSave = typename DenseVectorView<IndexType>::Engine(3, iSaveData, 1);
//
//  Decode and test the input parameters
//
    const bool wantBH = (job==TRSNA::Both);
    const bool wantS = (job==TRSNA::EigenvaluesOnly) || wantBH;
    const bool wantSP = (job==TRSNA::EigenvectorsOnly) || wantBH;

    const bool someCon = (howMany==TRSNA::Selected);
//
//  Set M to the number of eigenpairs for which condition numbers
//  are required, and test MM.
//
    if (someCon) {
        m = 0;
        bool pair = false;
        for (IndexType k=1; k<=n; ++k) {
            if (pair) {
                pair= false;
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
    } else {
        m = n;
    }

    if (mm<m) {
        ASSERT(0);
    }

    if (wantSP) {
        ASSERT(VL_.numCols()>=m);
        ASSERT(VR_.numCols()>=m);
    }
    // TODO: if one forgets to make this auto views const you get
    //       some error that is hard to understand for newbies ...
    //       Idea: disallow the creation of non-const views from
    //       const matrices/vectors.
    const auto VL = VL_(_,_(1,m));
    const auto VR = VR_(_,_(1,m));

//
//  Quick return if possible
//
    if (n==0) {
        return;
    }
//
    if (n==1) {
        if (someCon) {
            if (!select(1)) {
                return;
            }
        }
        if (wantS) {
            s(1) = One;
        }
        if (wantSP) {
            sep(1) = abs(T(1,1));
        }
        return;
    }
//
//  Get machine constants
//
    const ElementType eps = lamch<ElementType>(Precision);
    ElementType smallNum = lamch<ElementType>(SafeMin) / eps;
    ElementType bigNum = One/smallNum;
    labad(smallNum, bigNum);

    IndexType ks = 0;
    bool  pair = false;
    for (IndexType k=1; k<=n; ++k) {
//
//      Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
//
        if (pair) {
            pair = false;
            continue;
        } else {
            if (k<n) {
                pair = (T(k+1,k)!=Zero);
            }
        }
//
//      Determine whether condition numbers are required for the k-th
//      eigenpair.
//
        if (someCon) {
            if (pair) {
                if (!select(k) && !select(k+1)) {
                    continue;
                }
            } else {
                if (!select(k)) {
                    continue;
                }
            }
        }

        ++ks;

        if (wantS) {
//
//          Compute the reciprocal condition number of the k-th
//          eigenvalue.
//
            if (!pair) {
//
//              Real eigenvalue.
//
                const ElementType prod = VR(_,ks) * VL(_,ks);
                const ElementType rNrm = blas::nrm2(VR(_,ks));
                const ElementType lNrm = blas::nrm2(VL(_,ks));
                s(ks) = abs(prod) / (rNrm*lNrm);
            } else {
//
//              Complex eigenvalue.
//
                const ElementType prod1 = VR(_,ks) * VL(_,ks)
                                        + VR(_,ks+1) * VL(_,ks+1);
                const ElementType prod2 = VL(_,ks) * VR(_,ks+1)
                                        - VL(_,ks+1) * VR(_,ks);
                const ElementType rNrm = lapy2(blas::nrm2(VR(_,ks)),
                                               blas::nrm2(VR(_,ks+1)));
                const ElementType lNrm = lapy2(blas::nrm2(VL(_,ks)),
                                               blas::nrm2(VL(_,ks+1)));
                const ElementType cond = lapy2(prod1, prod2) / (rNrm*lNrm);
                s(ks)   = cond;
                s(ks+1) = cond;
            }
        }

        if (wantSP) {
//
//          Estimate the reciprocal condition number of the k-th
//          eigenvector.
//
//          Copy the matrix T to the array WORK and swap the diagonal
//          block beginning at T(k,k) to the (1,1) position.
//
            auto T_ = Work(_,_(1,n));
            T_ = T;
            IndexType iFirst =k;
            IndexType iLast = 1;
            IndexType iErr = trexc(false, T_, T_, iFirst, iLast, Work(_,n+1));

            ElementType est, mu, scale;
            IndexType   n2, nn;

            if (iErr==1 || iErr==2) {
//
//              Could not swap because blocks not well separated
//
                scale = One;
                est   = bigNum;
            } else {
//
//              Reordering successful
//
                if (Work(2,1)==Zero) {
//
//                  Form C = T22 - lambda*I in WORK(2:N,2:N).
//
                    for (IndexType i=2; i<=n; ++i) {
                        Work(i,i) -= Work(1,1);
                    }
                    n2 = 1;
                    nn = n - 1;
                } else {
//
//                  Triangularize the 2 by 2 block by unitary
//                  transformation U = [  cs   i*ss ]
//                                     [ i*ss   cs  ].
//                  such that the (1,1) position of WORK is complex
//                  eigenvalue lambda with positive imaginary part. (2,2)
//                  position of WORK is the complex eigenvalue lambda
//                  with negative imaginary  part.
//
                    mu = sqrt(abs(Work(1,2))) * sqrt(abs(Work(2,1)));
                    const ElementType delta = lapy2(mu, Work(2,1));
                    const ElementType cs = mu / delta;
                    const ElementType sn = -Work(2,1) / delta;
//
//                  Form
//
//                  C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
//                                           [   mu                     ]
//                                           [         ..               ]
//                                           [             ..           ]
//                                           [                  mu      ]
//                  where C**T is transpose of matrix C,
//                  and RWORK is stored starting in the N+1-st column of
//                  WORK.
//
                    Work(2,_(3,n)) *= cs;
                    for (IndexType j=3; j<=n; ++j) {
                        Work(j,j) -= Work(1,1);
                    }
                    Work(2,2) = Zero;

                    Work(1,n+1) =Two*mu;
                    for (IndexType i=2; i<=n-1; ++i) {
                        Work(i,n+1) = sn*Work(1,i+1);
                    }
                    n2 = 2;
                    nn = 2*(n-1);
                }
//
//              Estimate norm(inv(C**T))
//
                est = Zero;
                IndexType kase = 0;
                do {
                    auto v_ = Work(_,_(n+2,n+3)).vectorView(1,nn);
                    auto x_ = Work(_,_(n+4,n+5)).vectorView(1,nn);
                    auto iSgn_ = iWork(_(1,nn));
                    lacn2(v_, x_, iSgn_, est, kase, iSave);
                    if (kase==0) {
                        break;
                    } else {
                        auto T = Work(_(2,n),_(2,n));
                        auto b = Work(_(1,n-1),n+1);
                        auto x = Work(_,_(n+4,n+5)).vectorView(1,2*(n-1));
                        auto w = Work(_(1,n-1),n+6);
                        ElementType dummyMu{};
                        if (kase==1) {
                            if (n2==1) {
//
//                              Real eigenvalue: solve C**T*x = scale*c.
//
                                laqtr(true, true, T, b, dummyMu, scale, x, w);
                            } else {
//
//                              Complex eigenvalue: solve
//                              C**T*(p+iq) = scale*(c+id) in real arithmetic.
//
                                laqtr(true, false, T, b, mu, scale, x, w);
                            }
                        } else {
                            if (n2==1) {
//
//                              Real eigenvalue: solve C*x = scale*c.
//
                                laqtr(false, true, T, w, dummyMu, scale, x, w);
                            } else {
//
//                              Complex eigenvalue: solve
//                              C*(p+iq) = scale*(c+id) in real arithmetic.
//
                                laqtr(false, false, T, b, mu, scale, x, w);

                            }
                        }
                    }
                } while (true);
            }

            sep(ks) = scale / max(est, smallNum);
            if (pair) {
                sep(ks+1) = sep(ks);
            }
        }

        if (pair) {
            ++ks;
        }
    }
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename M, typename MM,
          typename MWORK, typename VIWORK>
void
trsna_impl(TRSNA::Job                    job,
           TRSNA::HowMany                howMany,
           const DenseVector<VSELECT>    &select,
           const GeMatrix<MT>            &T,
           const GeMatrix<MVL>           &VL,
           const GeMatrix<MVR>           &VR,
           DenseVector<VS>               &s,
           DenseVector<VSEP>             &sep,
           const MM                      &mm,
           M                             &m,
           GeMatrix<MWORK>               &Work,
           DenseVector<VIWORK>           &iWork)
{
    typedef typename GeMatrix<MT>::IndexType     IndexType;

    DenseVector<Array<IndexType> > select_ = select;

    cxxlapack::trsna<IndexType>(getF77Char(job),
                                getF77Char(howMany),
                                select_.data(),
                                T.numRows(),
                                T.data(),
                                T.leadingDimension(),
                                VL.data(),
                                VL.leadingDimension(),
                                VR.data(),
                                VR.leadingDimension(),
                                s.data(),
                                sep.data(),
                                mm,
                                m,
                                Work.data(),
                                Work.leadingDimension(),
                                iWork.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename MM, typename M,
          typename MWORK, typename VIWORK>
void
trsna(TRSNA::Job                    job,
      TRSNA::HowMany                howMany,
      const DenseVector<VSELECT>    &select,
      const GeMatrix<MT>            &T,
      const GeMatrix<MVL>           &VL,
      const GeMatrix<MVR>           &VR,
      DenseVector<VS>               &s,
      DenseVector<VSEP>             &sep,
      const MM                      &mm,
      M                             &m,
      GeMatrix<MWORK>               &Work,
      DenseVector<VIWORK>           &iWork)
{

#if defined(CHECK_CXXLAPACK) || !defined(NDEBUG)
    typedef typename GeMatrix<MT>::IndexType IndexType;
    const IndexType n = T.numRows();
#endif

#   ifndef NDEBUG
    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());

    if (howMany!=TRSNA::All) {
        ASSERT(select.firstIndex()==1);
        ASSERT(select.length()==n);
    }

    if (job!=TRSNA::EigenvectorsOnly) {
        ASSERT(VL.firstRow()==1);
        ASSERT(VL.firstCol()==1);
        ASSERT(VL.numRows()==n);
    }

    if (job!=TRSNA::EigenvectorsOnly) {
        ASSERT(VR.firstRow()==1);
        ASSERT(VR.firstCol()==1);
        ASSERT(VR.numRows()==n);
    }

    ASSERT(s.firstIndex()==1);
    ASSERT(s.length()==mm);

    ASSERT(sep.firstIndex()==1);
    ASSERT(sep.length()==mm);

    if (job!=TRSNA::EigenvaluesOnly) {
        ASSERT(Work.firstRow()==1);
        ASSERT(Work.firstCol()==1);
        ASSERT(Work.numRows()==n);
        ASSERT(Work.numCols()==n+6);

        ASSERT(iWork.firstIndex()==1);
        ASSERT(iWork.length()==2*(n-1));
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename DenseVector<VS>::NoView      s_org     = s;
    typename DenseVector<VSEP>::NoView    sep_org   = sep;
    M                                     m_org     = m;
    typename GeMatrix<MWORK>::NoView      Work_org  = Work;
    typename DenseVector<VIWORK>::NoView  iWork_org = iWork;
#   endif

    LAPACK_SELECT::trsna_impl(job, howMany, select, T, VL, VR,
                              s, sep, mm, m, Work, iWork);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename DenseVector<VS>::NoView      s_generic     = s;
    typename DenseVector<VSEP>::NoView    sep_generic   = sep;
    M                                     m_generic     = m;
    typename GeMatrix<MWORK>::NoView      Work_generic  = Work;
    typename DenseVector<VIWORK>::NoView  iWork_generic = iWork;

//
//  restore output arguments
//
    s       = s_org;
    sep     = sep_org;
    m       = m_org;
    Work    = Work_org;
    iWork   = iWork_org;

    external::trsna_impl(job, howMany, select, T, VL, VR,
                         s, sep, mm, m, Work, iWork);

    bool failed = false;
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

    if (! isIdentical(m_generic, m, "m_generic", "m")) {
        std::cerr << "CXXLAPACK: m_generic = " << m_generic << std::endl;
        std::cerr << "F77LAPACK: m = " << m << std::endl;
        failed = true;
    }

    if (! isIdentical(Work_generic, Work, "Work_generic", "Work")) {
        std::cerr << "CXXLAPACK: Work_generic = " << Work_generic << std::endl;
        std::cerr << "F77LAPACK: Work = " << Work << std::endl;
        failed = true;
    }

    if (! isIdentical(iWork_generic, iWork, "iWork_generic", "iWork")) {
        std::cerr << "CXXLAPACK: iWork_generic = "
                  << iWork_generic << std::endl;
        std::cerr << "F77LAPACK: iWork = " << iWork << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "n = " << n << std::endl;
        ASSERT(0);
    } else {
//        std::cerr << "passed: trsna.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename VSELECT, typename MT, typename MVL, typename MVR,
          typename VS, typename VSEP, typename MM, typename M,
          typename MWORK, typename VIWORK>
void
trsna(TRSNA::Job                    job,
      TRSNA::HowMany                howMany,
      const VSELECT                 &select,
      const MT                      &T,
      const MVL                     &VL,
      const MVR                     &VR,
      VS                            &&s,
      VSEP                          &&sep,
      const MM                      &mm,
      M                             &&m,
      MWORK                         &&Work,
      VIWORK                        &&iWork)
{
    CHECKPOINT_ENTER;
    trsna(job, howMany, select, T, VL, VR, s, sep, mm, m, Work, iWork);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_TRSNA_TCC
