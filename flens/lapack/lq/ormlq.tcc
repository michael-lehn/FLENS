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
       SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LQ_ORMLQ_TCC
#define FLENS_LAPACK_LQ_ORMLQ_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MA, typename MC>
typename GeMatrix<MC>::IndexType
ormlq_generic_wsq(Side              side,
                  Transpose         trans,
                  GeMatrix<MA>      &A,
                  GeMatrix<MC>      &C)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    typedef typename GeMatrix<MC>::View         GeView;
    typedef typename GeView::Engine             GeViewEngine;

//
//  Paramter for maximum block size and buffer for TrMatrix Tr.
//
    const IndexType nbMax   = 64;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    const IndexType nw      = (side==Left) ? n : m;

//
//  Determine the block size.  nb may be at most nbMax, where nbMax
//  is used to define the local array tr.
//
    char opt[3];
    opt[0] = char(side);
    if (trans==NoTrans) {
        opt[1] = 'N';
    } else if (trans==Conj) {
        opt[1] = 'R';
    } else if (trans==Trans) {
        opt[1] = 'T';
    } else if (trans==ConjTrans) {
        opt[1] = 'C';
    }
    opt[2] = 0;

    IndexType nb = min(nbMax, IndexType(ilaenv<T>(1, "ORMLQ", opt, m, n, k)));
    return max(IndexType(1), nw)*nb;
}

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormlq_generic(Side                      side,
              Transpose                 trans,
              GeMatrix<MA>              &A,
              const DenseVector<VTAU>   &tau,
              GeMatrix<MC>              &C,
              DenseVector<VWORK>        &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MC>::ElementType  T;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    typedef typename GeMatrix<MC>::View         GeView;
    typedef typename GeView::Engine             GeViewEngine;

    const Underscore<IndexType> _;


//
//  Paramter for maximum block size and buffer for TrMatrix Tr.
//
    const IndexType     nbMax = 64;
    const IndexType     ldt = nbMax + 1;
    T                   trBuffer[nbMax*ldt];

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    const bool noTrans = (trans==Trans || trans==ConjTrans) ? false : true;
//
//  nq is the order of Q and nw is the minimum dimension of work
//
    IndexType nq, nw;
    if (side==Left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
//
//  Determine the block size.  nb may be at most nbMax, where nbMax
//  is used to define the local array tr.
//
    char opt[3];
    opt[0] = (side==Left) ? 'L' : 'R';
    if (trans==NoTrans) {
        opt[1] = 'N';
    } else if (trans==Conj) {
        opt[1] = 'R';
    } else if (trans==Trans) {
        opt[1] = 'T';
    } else if (trans==ConjTrans) {
        opt[1] = 'C';
    }
    opt[2] = 0;

    IndexType nb = min(nbMax, IndexType(ilaenv<T>(1, "ORMLQ", opt, m, n, k)));
    IndexType lWorkOpt = max(IndexType(1), nw)*nb;

    if (work.length()==0) {
        work.resize(lWorkOpt);
    }

//
//  Quick return if possible
//
    if ((m==0) || (n==0) || (k==0)) {
        work(1) = 1;
        return;
    }

    IndexType nbMin = 2;
    IndexType iws;
    if ((nb>1) && (nb<k)) {
        iws = lWorkOpt;
        if (work.length()<iws) {
            nb = work.length()/nw;
            nbMin = max(nbMin, IndexType(ilaenv<T>(2, "ORMLQ", opt, m, n, k)));
        }
    } else {
        iws = nw;
    }

    if ((nb<nbMin) || (nb>=k)) {
//
//      Use unblocked code
//
        auto _work = (side==Left) ? work(_(1,n)) : work(_(1,m));
        orml2(side, trans, A, tau, C, _work);
    } else {
//
//      Use blocked code
//
        IndexType iBeg, iInc, iEnd;
        if ((side==Left && noTrans) || (side==Right && !noTrans)) {
            iBeg = 1;
            iEnd = ((k-1)/nb)*nb + 1;
            iInc = nb;
        } else {
            iBeg = ((k-1)/nb)*nb + 1;
            iEnd = 1;
            iInc = -nb;
        }
        iEnd += iInc;

        IndexType ic, jc;
        if (side==Left) {
            jc = 1;
        } else {
            ic = 1;
        }

        const Transpose transT = (trans==NoTrans) ? Trans : NoTrans;

        typename GeMatrix<MA>::View Work(nw, nb, work);

        for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
            const IndexType ib = min(nb, k-i+1);
            GeView          Tr = GeViewEngine(ib, ib, trBuffer, ldt);
//
//          Form the triangular factor of the block reflector
//          H = H(i) H(i+1) . . . H(i+ib-1)
//
            larft(Forward, RowWise, nq-i+1,
                  A(_(i,i+ib-1),_(i,nq)), tau(_(i,i+ib-1)), Tr.upper());

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

            larfb(side, transT,
                  Forward, RowWise,
                  A(_(i,i+ib-1),_(i,nq)), Tr.upper(),
                  C(_(ic,m),_(jc,n)),
                  Work(_,_(1,ib)));
        }
    }
    work(1) = lWorkOpt;
}

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename MC>
typename GeMatrix<MC>::IndexType
ormlq_wsq(Side              side,
          Transpose         trans,
          GeMatrix<MA>      &A,
          GeMatrix<MC>      &C)
{
    typedef typename GeMatrix<MC>::IndexType    IndexType;
    typedef typename GeMatrix<MC>::ElementType  T;

    T                   WORK, DUMMY;
    const IndexType     LWORK   = -1;

    cxxlapack::ormlq<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                A.numRows(),
                                A.data(),
                                A.leadingDimension(),
                                &DUMMY,
                                C.data(),
                                C.leadingDimension(),
                                &WORK,
                                LWORK);
    return WORK;
}

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormlq(Side                       side,
      Transpose                  trans,
      GeMatrix<MA>               &A,
      const DenseVector<VTAU>    &tau,
      GeMatrix<MC>               &C,
      DenseVector<VWORK>         &work)
{
    typedef typename GeMatrix<MC>::IndexType  IndexType;

    cxxlapack::ormlq<IndexType>(getF77Char(side),
                                getF77Char(trans),
                                C.numRows(),
                                C.numCols(),
                                A.numRows(),
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

template <typename MA, typename MC>
typename GeMatrix<MC>::IndexType
ormlq_wsq(Side              side,
          Transpose         trans,
          GeMatrix<MA>      &A,
          GeMatrix<MC>      &C)
{
    typedef typename GeMatrix<MC>::IndexType    IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    if (side==Left) {
        ASSERT(A.numCols()==m);
    } else {
        ASSERT(A.numCols()==n);
    }
#   endif

//
//  Call implementation
//
    const IndexType info = ormlq_generic(side, trans, A, C);

#   ifdef CHECK_CXXLAPACK
//
//  Compare generic results with results from the native implementation
//
    const IndexType _info = ormlq_native(side, trans, A, C);

    ASSERT(info==_info);
#   endif
    return info;
}

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormlq(Side                      side,
      Transpose                 trans,
      GeMatrix<MA>              &A,
      const DenseVector<VTAU>   &tau,
      GeMatrix<MC>              &C,
      DenseVector<VWORK>        &work)
{
//
//  Test the input parameters
//
#   ifndef NDEBUG
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = A.numRows();

    ASSERT(tau.length()==k);

    if (side==Left) {
        ASSERT(A.numCols()==m);
    } else {
        ASSERT(A.numCols()==n);
    }

    if (work.length()>0) {
        if (side==Left) {
            ASSERT(work.length()>=n);
        } else {
            ASSERT(work.length()>=m);
        }
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org      = A;
    typename GeMatrix<MC>::NoView       C_org      = C;
    typename DenseVector<VWORK>::NoView work_org   = work;
#   endif

//
//  Call implementation
//
    ormlq_generic(side, trans, A, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic       = A;
    typename GeMatrix<MC>::NoView       C_generic       = C;
    typename DenseVector<VWORK>::NoView work_generic    = work;

//
//  restore output arguments
//
    A = A_org;
    C = C_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    ormlq_native(side, trans, A, tau, C, work);

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
        std::cerr << "error in: ormlq.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: ormlq.tcc" << std::endl;
    }
#   endif
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename MC>
typename MC::IndexType
ormlq_wsq(Side          side,
          Transpose     trans,
          MA            &&A,
          MC            &&C)
{
    typedef typename MC::IndexType IndexType;

    CHECKPOINT_ENTER;
    const IndexType info = ormlq_wsq(side, trans, A, C);
    CHECKPOINT_LEAVE;

    return info;
}

template <typename MA, typename VTAU, typename MC, typename VWORK>
void
ormlq(Side              side,
      Transpose         trans,
      MA                &&A,
      const VTAU        &tau,
      MC                &&C,
      VWORK             &&work)
{
    // TODO: asser that A is non-const
    CHECKPOINT_ENTER;
    ormlq(side, trans, A, tau, C, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LQ_ORMLQ_TCC
