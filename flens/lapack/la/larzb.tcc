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
       SUBROUTINE DLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,
      $                   LDV, T, LDT, C, LDC, WORK, LDWORK )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LARZB_TCC
#define FLENS_LAPACK_LA_LARZB_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- larzb [real variant] ------------------------------------------------------

template <typename MV, typename MT, typename MC, typename MWORK>
void
larzb_impl(Side                     side,
           Transpose                trans,
           Direction                direction,
           StoreVectors             storeVectors,
           const GeMatrix<MV>       &V,
           const TrMatrix<MT>       &T,
           GeMatrix<MC>             &C,
           GeMatrix<MWORK>          &Work)
{
    typedef typename GeMatrix<MC>::ElementType  ElementType;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const ElementType  One(1);

    const Underscore<IndexType>  _;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = T.dim();
    const IndexType l = (storeVectors==RowWise) ? V.numCols() : V.numRows();

//
//  Quick return if possible
//
    if (m==0 || n==0) {
        return;
    }

//
//  Check for currently supported options
//
    ASSERT(direction==Backward);
    ASSERT(storeVectors==RowWise);

    Transpose transT = (trans==NoTrans) ? Trans : NoTrans;

    if (side==Left) {
//
//      Form  H * C  or  H**T * C
//
//      W( 1:n, 1:k ) = C( 1:k, 1:n )**T
//
        for (IndexType j=1; j<=k; ++j) {
            Work(_,j) = C(j,_);
        }
//
//      W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
//                      C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T
//
        if (l>0) {
            blas::mm(Trans, Trans, One, C(_(m-l+1,m),_), V, One, Work);
        }

//
//      W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:n, 1:k ) * T
//
        blas::mm(Right, transT, One, T, Work);

//
//      C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T
//
        for (IndexType j=1; j<=n; ++j) {
            for (IndexType i=1; i<=k; ++i) {
                C(i,j) -= Work(j,i);
            }
        }

//
//      C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
//                          V( 1:k, 1:l )**T * W( 1:n, 1:k )**T
//
        if (l>0) {
            blas::mm(Trans, Trans, -One, V, Work, One, C(_(m-l+1,m),_));
        }

    } else if (side==Right) {
//
//      Form  C * H  or  C * H**T
//
//      W( 1:m, 1:k ) = C( 1:m, 1:k )
//
        Work = C(_,_(1,k));

//
//      W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
//                      C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T
//
        if (l>0) {
            blas::mm(NoTrans, Trans, One, C(_,_(n-l+1,n)), V, One, Work);
        }

//
//      W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T
//
        blas::mm(Right, trans, One, T, Work);

//
//      C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
//
        C(_,_(1,k)) -= Work;

//
//      C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
//                          W( 1:m, 1:k ) * V( 1:k, 1:l )
//
        if (l>0) {
            blas::mm(NoTrans, NoTrans, -One, Work, V, One, C(_,_(n-l+1,n)));
        }

    }

}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- larzb [real variant] ------------------------------------------------------

template <typename MV, typename MT, typename MC, typename MWORK>
void
larzb_impl(Side                     side,
           Transpose                trans,
           Direction                direction,
           StoreVectors             storeVectors,
           const GeMatrix<MV>       &V,
           const TrMatrix<MT>       &T,
           GeMatrix<MC>             &C,
           GeMatrix<MWORK>          &Work)
{
    cxxlapack::larzb(getF77Char(side),
                     getF77Char(trans),
                     getF77Char(direction),
                     getF77Char(storeVectors),
                     C.numRows(),
                     C.numCols(),
                     T.dim(),
                     (storeVectors==RowWise) ? V.numCols() : V.numRows(),
                     V.data(),
                     V.leadingDimension(),
                     T.data(),
                     T.leadingDimension(),
                     C.data(),
                     C.leadingDimension(),
                     Work.data(),
                     Work.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- larzb[real variant] ------------------------------------------------------

template <typename MV, typename MT, typename MC, typename MWORK>
typename RestrictTo<IsRealGeMatrix<MV>::value
                 && IsRealTrMatrix<MT>::value
                 && IsRealGeMatrix<MC>::value
                 && IsRealGeMatrix<MWORK>::value,
         void>::Type
larzb(Side           side,
      Transpose      trans,
      Direction      direction,
      StoreVectors   storeVectors,
      const MV       &V,
      const MT       &T,
      MC             &&C,
      MWORK          &&Work)
{

//
//  Test the input parameters
//
#   ifndef NDEBUG

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename MatrixC::IndexType     IndexType;
    
//  From LAPACK:  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
    ASSERT(storeVectors==RowWise);
    ASSERT(direction==Backward);

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = T.dim();
    const IndexType l = (storeVectors==RowWise) ? V.numCols() : V.numRows();

    if (side==Left) {
        ASSERT(m>=l);
        ASSERT(Work.numRows()==n);
    } else {
        ASSERT(n>=l);
        ASSERT(Work.numRows()==m);
    }

    if (storeVectors==RowWise) {
        ASSERT(V.numRows()==k);
    } else {
        ASSERT(V.numCols()==k);
    }
    ASSERT(T.upLo()==Lower);
    ASSERT(Work.numCols()==k);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<MWORK>::Type MatrixWork;
    
    typename MatrixC::NoView        C_org    = C;
    typename MatrixWork::NoView     Work_org = Work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::larzb_impl(side, trans, direction, storeVectors,
                              V, T, C, Work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixC::NoView        C_generic    = C;
    typename MatrixWork::NoView     Work_generic = Work;

//
//  restore output arguments
//
    C    = C_org;
    Work = Work_org;

//
//  Compare generic results with results from the native implementation
//
    external::larzb_impl(side, trans, direction, storeVectors, V, T, C, Work);

    bool failed = false;
    if (! isIdentical(C_generic, C, "C_generic", "C")) {
        std::cerr << "CXXLAPACK: C_generic = " << C_generic << std::endl;
        std::cerr << "F77LAPACK: C = " << C << std::endl;
        failed = true;
    }
    if (! isIdentical(Work_generic, Work, "Work_generic", "Work")) {
        std::cerr << "CXXLAPACK: Work_generic = " << Work_generic << std::endl;
        std::cerr << "F77LAPACK: Work = " << Work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: larzb.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: larzb.tcc" << std::endl;
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LARZB_TCC
