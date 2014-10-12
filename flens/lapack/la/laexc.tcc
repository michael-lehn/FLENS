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
       SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
      $                   INFO )
 *
 *  -- LAPACK auxiliary routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 */

#ifndef FLENS_LAPACK_LA_LAEXC_TCC
#define FLENS_LAPACK_LA_LAEXC_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
laexc_impl(bool                          computeQ,
           GeMatrix<MT>                  &T,
           GeMatrix<MQ>                  &Q,
           IndexType                     j1,
           IndexType                     n1,
           IndexType                     n2,
           DenseVector<VWORK>            &work)
{
    using std::abs;
    using flens::max;

    typedef typename GeMatrix<MT>::View         GeMatrixView;
    typedef typename GeMatrix<MT>::VectorView   DenseVectorView;
    typedef typename GeMatrix<MT>::ElementType  ElementType;

    const ElementType Zero(0), One(1), Ten(10);

    const IndexType n = T.numRows();

    const Underscore<IndexType>  _;

//
//    .. Local Arrays ..
//
    ElementType dData_[16], xData_[4];
    GeMatrixView  D = typename GeMatrixView::Engine(4, 4, dData_, 4);
    GeMatrixView  X = typename GeMatrixView::Engine(2, 2, xData_, 2);

    ElementType uData_[3], u1Data_[3], u2Data_[3];
    DenseVectorView u  = typename DenseVectorView::Engine(3, uData_);
    DenseVectorView u1 = typename DenseVectorView::Engine(3, u1Data_);
    DenseVectorView u2 = typename DenseVectorView::Engine(3, u2Data_);
//
//  Quick return if possible
//
    if (n==0 || n1==0 || n2==0) {
        return 0;
    }
    if (j1+n1>n) {
        return 0;
    }

    const IndexType j2 = j1 + 1;
    const IndexType j3 = j1 + 2;
    const IndexType j4 = j1 + 3;

    ElementType t11, t22, t33;

    if (n1==1 && n2==1) {
//
//      Swap two 1-by-1 blocks.
//
        t11 = T(j1,j1);
        t22 = T(j2,j2);
//
//      Determine the transformation to perform the interchange.
//
        ElementType cs, sn, temp;
        lartg(T(j1,j2), t22-t11, cs, sn, temp);
//
//      Apply transformation to the matrix T.
//
        if (j3<=n) {
            blas::rot(T(j1,_(j3,n)), T(j2,_(j3,n)), cs, sn);
        }
        blas::rot(T(_(1,j1-1),j1), T(_(1,j1-1),j2), cs, sn);

        T(j1,j1) = t22;
        T(j2,j2) = t11;

        if (computeQ) {
//
//          Accumulate transformation in the matrix Q.
//
            blas::rot(Q(_,j1), Q(_,j2), cs, sn);
        }

    } else {
//
//      Swapping involves at least one 2-by-2 block.
//
//      Copy the diagonal block of order N1+N2 to the local array D
//      and compute its norm.
//
        const IndexType nd = n1 + n2;
        auto D_ = D(_(1,nd),_(1,nd));

        D_ = T(_(j1,j1+nd-1),_(j1,j1+nd-1));
        ElementType normD = lan(MaximumNorm, D_);

        ElementType cs, sn, wr1, wr2, wi1, wi2;

        ElementType scale, normX, tau, tau1, tau2;
//
//      Compute machine-dependent threshold for test for accepting
//      swap.
//
        const ElementType eps = lamch<ElementType>(Precision);
        const ElementType smallNum = lamch<ElementType>(SafeMin) / eps;
        const ElementType thresh = max(Ten*eps*normD, smallNum);
//
//      Solve T11*X - X*T22 = scale*T12 for X.
//
        const auto T11 = D(_(1,n1),_(1,n1));
        const auto T12 = D(_(1,n1),_(n1+1,nd));
        const auto T22 = D(_(n1+1,nd),_(n1+1,nd));

        auto X_ = X(_(1,n1),_(1,n2));

        lasy2(false, false, IndexType(-1), T11, T22, T12, scale, X_, normX);
//
//      Swap the adjacent diagonal blocks.
//
        const IndexType k = n1 + n1 + n2 - 3;

        switch (k) {
//
//          N1 = 1, N2 = 2: generate elementary reflector H so that:
//
//          ( scale, X11, X12 ) H = ( 0, 0, * )
//
        case 1:
            u(1) = scale;
            u(2) = X(1,1);
            u(3) = X(1,2);
            larfg(IndexType(3), u(3), u(_(1,2)), tau);
            u(3) = One;
            t11 = T(j1,j1);
//
//          Perform swap provisionally on diagonal block in D.
//
            larfx(Left, u, tau, D_, work(_(1,3)));
            larfx(Right, u, tau, D_, work(_(1,3)));
//
//          Test whether to reject swap.
//
            if (max(abs(D(3,1)), abs(D(3,2)), abs(D(3,3)-t11))>thresh) {
//
//              Return 1 if swap was rejected.
//
                return 1;
            }
//
//          Accept swap: apply transformation to the entire matrix T.
//
            larfx(Left, u, tau, T(_(j1,j1+3-1),_(j1,n)), work(_(1,n-j1+1)));
            larfx(Right, u, tau, T(_(1,j2),_(j1,j1+3-1)), work(_(1,j2)));

            T(j3,j1) = Zero;
            T(j3,j2) = Zero;
            T(j3,j3) = t11;

            if (computeQ) {
//
//              Accumulate transformation in the matrix Q.
//
                larfx(Right, u, tau, Q(_,_(j1,j1+3-1)), work);
            }
            break;

        case 2:
//
//          N1 = 2, N2 = 1: generate elementary reflector H so that:
//
//          H (  -X11 ) = ( * )
//            (  -X21 ) = ( 0 )
//            ( scale ) = ( 0 )
//
            u(1) = -X(1,1);
            u(2) = -X(2,1);
            u(3) = scale;
            larfg(IndexType(3), u(1), u(_(2,3)), tau);
            u(1) = One;
            t33 = T(j3,j3);
//
//          Perform swap provisionally on diagonal block in D.
//
            larfx(Left, u, tau, D(_(1,3),_(1,3)), work(_(1,3)));
            larfx(Right, u, tau, D(_(1,3),_(1,3)), work(_(1,3)));
//
//          Test whether to reject swap.
//
            if (max(abs(D(2,1)), abs(D(3,1)), abs(D(1,1)-t33))>thresh) {
//
//              Return 1 if swap was rejected.
//
                return 1;
            }
//
//          Accept swap: apply transformation to the entire matrix T.
//
            larfx(Right, u, tau, T(_(1,j3),_(j1, j1+3-1)), work(_(1,j3)));
            larfx(Left, u, tau, T(_(j1,j1+3-1),_(j2,n)), work(_(1,n-j1)));

            T(j1,j1) = t33;
            T(j2,j1) = Zero;
            T(j3,j1) = Zero;

            if (computeQ) {
//
//              Accumulate transformation in the matrix Q.
//
                larfx(Right, u, tau, Q(_,_(j1,j1+3-1)), work);
            }
            break;

        case 3:
//
//          N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
//          that:
//
//          H(2) H(1) (  -X11  -X12 ) = (  *  * )
//                    (  -X21  -X22 )   (  0  * )
//                    ( scale    0  )   (  0  0 )
//                    (    0  scale )   (  0  0 )
//
            u1(1) = -X(1,1);
            u1(2) = -X(2,1);
            u1(3) = scale;
            larfg(IndexType(3), u1(1), u1(_(2,3)), tau1);
            u1(1) = One;

            const ElementType temp = -tau1*(X(1,2)+u1(2)*X(2,2));
            u2(1) = -temp*u1(2) - X(2,2);
            u2(2) = -temp*u1(3);
            u2(3) = scale;
            larfg(IndexType(3), u2(1), u2(_(2,3)), tau2);
            u2(1) = One;
//
//          Perform swap provisionally on diagonal block in D.
//
            larfx(Left,  u1, tau1, D(_(1,3),_(1,4)), work(_(1,4)));
            larfx(Right, u1, tau1, D(_(1,4),_(1,3)), work(_(1,4)));
            larfx(Left,  u2, tau2, D(_(2,4),_(1,4)), work(_(1,4)));
            larfx(Right, u2, tau2, D(_(1,4),_(2,4)), work(_(1,4)));
//
//          Test whether to reject swap.
//
            if (max(abs(D(3,1)), abs(D(3,2)), abs(D(4,1)), abs(D(4,2)))>thresh)
            {
//
//              Return 1 if swap was rejected.
//
                return 1;
            }
//
//          Accept swap: apply transformation to the entire matrix T.
//
            larfx(Left,  u1, tau1, T(_(j1,j1+3-1),_(j1,n)), work(_(1,n-j1+1)));
            larfx(Right, u1, tau1, T(_(1,j4),_(j1,j1+3-1)), work(_(1,j4)));
            larfx(Left,  u2, tau2, T(_(j2,j2+3-1),_(j1,n)), work(_(1,n-j1+1)));
            larfx(Right, u2, tau2, T(_(1,j4),_(j2,j2+3-1)), work(_(1,j4)));

            T(j3,j1) = Zero;
            T(j3,j2) = Zero;
            T(j4,j1) = Zero;
            T(j4,j2) = Zero;

            if (computeQ) {
//
//              Accumulate transformation in the matrix Q.
//
                larfx(Right, u1, tau1, Q(_,_(j1,j1+3-1)), work);
                larfx(Right, u2, tau2, Q(_,_(j2,j2+3-1)), work);
            }
        }

        if (n2==2) {
//
//          Standardize new 2-by-2 block T11
//
            lanv2(T(j1,j1), T(j1,j2), T(j2,j1), T(j2,j2),
                  wr1, wi1, wr2, wi2, cs, sn);
            blas::rot(T(j1,_(j1+2,n)), T(j2,_(j1+2,n)), cs, sn);
            blas::rot(T(_(1,j1-1),j1), T(_(1,j1-1),j2), cs, sn);
            if (computeQ) {
                blas::rot(Q(_,j1), Q(_,j2), cs, sn);
            }
        }

        if (n1==2) {
//
//          Standardize new 2-by-2 block T22
//
            const IndexType j3 = j1 + n2;
            const IndexType j4 = j3 + 1;
            lanv2(T(j3,j3), T(j3,j4), T(j4,j3), T(j4,j4),
                  wr1, wi1, wr2, wi2, cs, sn);
            if (j3+2<=n) {
                blas::rot(T(j3,_(j3+2,n)), T(j4,_(j3+2,n)), cs, sn);
            }
            blas::rot(T(_(1,j3-1),j3), T(_(1,j3-1),j4), cs, sn);
            if (computeQ) {
                blas::rot(Q(_,j3), Q(_,j4), cs, sn);
            }
        }
    }
    return 0;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
laexc_impl(bool                          computeQ,
           GeMatrix<MT>                  &T,
           GeMatrix<MQ>                  &Q,
           IndexType                     j1,
           IndexType                     n1,
           IndexType                     n2,
           DenseVector<VWORK>            &work)
{
    IndexType  info;
    info = cxxlapack::laexc<IndexType>(computeQ,
                                       T.numRows(),
                                       T.data(),
                                       T.leadingDimension(),
                                       Q.data(),
                                       Q.leadingDimension(),
                                       j1,
                                       n1,
                                       n2,
                                       work.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
laexc(bool                          computeQ,
      GeMatrix<MT>                  &T,
      GeMatrix<MQ>                  &Q,
      IndexType                     j1,
      IndexType                     n1,
      IndexType                     n2,
      DenseVector<VWORK>            &work)
{
    LAPACK_DEBUG_OUT("BEGIN: laexc");

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
        ASSERT(Q.numRows()==Q.numCols());
        ASSERT(Q.numRows()==n);
    }

    ASSERT(j1>=1);
    ASSERT((n1==0) || (n1==1) || (n1==2));
    ASSERT((n2==0) || (n2==1) || (n2==2));

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MT>::NoView           T_org      = T;
    typename GeMatrix<MQ>::NoView           Q_org      = Q;
    typename DenseVector<VWORK>::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::laexc_impl(computeQ, T, Q,
                                               j1, n1, n2, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MT>::NoView           T_generic       = T;
    typename GeMatrix<MQ>::NoView           Q_generic       = Q;
    typename DenseVector<VWORK>::NoView     work_generic    = work;

//
//  restore output arguments
//
    T = T_org;
    Q = Q_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//

    IndexType info_ = external::laexc_impl(computeQ, T, Q, j1, n1, n2, work);

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

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(info, info_, " info", "info_")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: info_ = " << info_ << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    LAPACK_DEBUG_OUT("END: laexc");

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
laexc(bool              computeQ,
      MT                &&T,
      MQ                &&Q,
      IndexType         j1,
      IndexType         n1,
      IndexType         n2,
      VWORK             &&work)
{
    CHECKPOINT_ENTER;
    const IndexType info = laexc(computeQ, T, Q, j1, n1, n2, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAEXC_TCC
