/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_ORMRZ_TCC
#define FLENS_LAPACK_IMPL_ORMRZ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- ormrz [real variant] ------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
void
ormrz_impl(Side                     side,
           Transpose                trans,
           IndexType                l,
           GeMatrix<MA>             &A,
           const DenseVector<VTAU>  &tau,
           GeMatrix<MC>             &C,
           DenseVector<VWORK>       &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MC>::ElementType  ElementType;
    typedef typename GeMatrix<MC>::View         GeView;
    typedef typename GeView::Engine             GeViewEngine;

//
//  Paramter for maximum block size and buffer for TrMatrix Tr.
//
    const IndexType nbMax = 64;
    const IndexType ldTr  = nbMax + 1;
    ElementType     trBuffer[nbMax*ldTr];

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = tau.length();

    const Underscore<IndexType>  _;

    const IndexType nw = (side==Left) ? max(1, n) : max(1, m);

    IndexType nb, lWorkOpt;
    if (m==0 || n==0) {
        nb       = 1;
        lWorkOpt = 1;
    } else {
//
//      Determine the block size.  NB may be at most NBMAX, where
//      NBMAX is used to define the local array T.
//
        char opt[3];
        opt[0] = getF77Char(side);
        opt[1] = getF77Char(trans);
        opt[2] = 0;
        nb = min(nbMax, ilaenv<ElementType>(1, "ORMRQ", opt, m, n, k));
        lWorkOpt = nw*nb;
    }

    if (work.length()==0) {
        work.resize(lWorkOpt);
    }
    work(1) = lWorkOpt;

//
//  Quick return if possible
//
    if (m==0 || n==0) {
        work(1) = 1;
        return;
    }

    IndexType nbMin  = 2;
    IndexType ldWork = nw;
    if (nb>1 && nb<k) {
        if (work.length()<nw*nb) {
            char opt[3];
            opt[0] = getF77Char(side);
            opt[1] = getF77Char(trans);
            opt[2] = 0;

            nb = work.length() / ldWork;
            nbMin = min(nbMin, ilaenv<ElementType>(2, "ORMRQ", opt, m, n, k));
        }
    }

    if (nb<nbMin || nb>=k) {
//
//      Use unblocked code
//
        auto work_ = (side==Left) ? work(_(1,n)) : work(_(1,m));
        ormr3(side, trans, l, A, tau, C, work_);
    } else {
//
//      Use blocked code
//
        IndexType iBeg, iInc, iEnd;
        if ((side==Left && trans!=NoTrans)
         || (side==Right && trans==NoTrans)) {
            iBeg = 1;
            iEnd = ((k-1)/nb)*nb + 1;
            iInc = nb;
        } else {
            iBeg = ((k-1)/nb)*nb + 1;
            iEnd = 1;
            iInc = -nb;
        }
        iEnd += iInc;

        IndexType ja, ic, jc;
        if (side==Left) {
            ja = m-l+1;
            jc = 1;
        } else {
            ja = n-l+1;
            ic = 1;
        }

        typename GeMatrix<MA>::View Work(nw, nb, work, ldWork);

        Transpose  transT = (trans==NoTrans) ? Trans : NoTrans;

        for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
            const IndexType ib = min(nb, k-i+1);
            GeView          Tr = GeViewEngine(ib, ib, ldTr, trBuffer);
//
//          Form the triangular factor of the block reflector
//          H = H(i+ib-1) . . . H(i+1) H(i)
//
            auto       A_   = A(_(i,i+ib-1),_(ja,ja+l-1));
            const auto tau_ = tau(_(i,i+ib-1));
            larzt(Backward, RowWise, A_, tau_, Tr.lower());

            if (side==Left) {
//
//              H or H**T is applied to C(i:m,1:n)
//
                ic = i;
            } else {
//
//              H or H**T is applied to C(1:m,i:n)
//
                jc = i;
            }
//
//          Apply H or H**T
//
            auto       C_    = C(_(ic,m),_(jc,n));
            auto       Work_ = Work(_,_(1,ib));
            larzb(side, transT, Backward, RowWise, A_, Tr.lower(), C_, Work_);
        }

    }

    work(1) = lWorkOpt;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- ormrz [real variant] ------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
void
ormrz_impl(Side                     side,
           Transpose                trans,
           IndexType                l,
           GeMatrix<MA>             &A,
           const DenseVector<VTAU>  &tau,
           GeMatrix<MC>             &C,
           DenseVector<VWORK>       &work)
{
    cxxlapack::ormrz<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                tau.length(),
                                l,
                                A.data(),
                                A.leadingDimension(),
                                tau.data(),
                                C.data(),
                                C.leadingDimension(),
                                work.data(),
                                work.length());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- ormrz [real variant] ------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
typename RestrictTo<IsInteger<IndexType>::value
                 && IsRealGeMatrix<MA>::value
                 && IsRealDenseVector<VTAU>::value
                 && IsRealGeMatrix<MC>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
ormrz(Side          side,
      Transpose     trans,
      IndexType     l,
      MA            &&A,
      const VTAU    &tau,
      MC            &&C,
      VWORK         &&work)
{
    using std::max;
    using std::min;

    LAPACK_DEBUG_OUT("ormrz [real]");

//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename RemoveRef<VWORK>::Type VectorWork;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = tau.length();

    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    ASSERT(l>=0);
    if (side==Left) {
        ASSERT(m>=k);
        ASSERT(m>=l);
    } else {
        ASSERT(n>=k);
        ASSERT(n>=l);
    }

    ASSERT(A.numRows()==k);
    if (side==Left) {
        ASSERT(A.numCols()==m);
    } else {
        ASSERT(A.numCols()==n);
    }

    ASSERT(tau.firstIndex()==1);

    ASSERT(C.firstRow()==1);
    ASSERT(C.firstCol()==1);

    ASSERT(work.firstIndex()==1);

    if (work.length()>0) {
        if (side==Left) {
            ASSERT(work.length()>=max(IndexType(1),n));
        } else {
            ASSERT(work.length()>=max(IndexType(1),m));
        }
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixA::NoView        A_org      = A;
    typename MatrixC::NoView        C_org      = C;
    typename VectorWork::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::ormrz_impl(side, trans, l, A, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView        A_generic    = A;
    typename MatrixC::NoView        C_generic    = C;
    typename VectorWork::NoView     work_generic = work;

//
//  restore output arguments
//
    A    = A_org;
    C    = C_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::ormrz_impl(side, trans, l, A, tau, C, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
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
        std::cerr << "error in: ormrz.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: ormrz.tcc" << std::endl;
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_ORMRZ_TCC
