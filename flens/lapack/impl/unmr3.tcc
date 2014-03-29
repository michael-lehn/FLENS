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
       SUBROUTINE ZUNMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
      $                   WORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_IMPL_UNMR3_TCC
#define FLENS_LAPACK_IMPL_UNMR3_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- unmr3 [complex variant] ---------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
void
unmr3_impl(Side                      side,
           Transpose                 trans,
           IndexType                 l,
           const GeMatrix<MA>        &A,
           const DenseVector<VTAU>   &tau,
           GeMatrix<MC>              &C,
           DenseVector<VWORK>        &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;

    const IndexType k =  tau.length();
    const IndexType m =  C.numRows();
    const IndexType n =  C.numCols();

    const IndexType na = A.numCols();
    const IndexType ja = na-l+1;

    ElementType taui(0);

    const Underscore<IndexType>  _;

//
//  Quick return if possible
//
    if (m==0 || n==0 || k==0) {
        return;
    }

    IndexType iBeg, iInc, iEnd;
    if ((side==Left && trans!=NoTrans) || (side==Right && trans==NoTrans)) {
        iBeg = 1;
        iEnd = k;
        iInc = 1;
    } else {
        iBeg = k;
        iEnd = 1;
        iInc = -1;
    }
    iEnd += iInc;

    IndexType ic, jc;
    if (side==Left) {
        jc = 1;
    } else {
        ic = 1;
    }

    for (IndexType i=iBeg; i!=iEnd; i+=iInc) {
        if (side==Left) {
//
//          H(i) or H(i)**H is applied to C(i:m,1:n)
//
            ic = i;
        } else {
//
//          H(i) or H(i)**H is applied to C(1:m,i:n)
//
            jc = i;
        }
//
//      Apply H(i) or H(i)**H
//
        if (trans==NoTrans) {
            taui = tau(i);
        } else {
            taui = conj(tau(i));
        }
        larz(side, A(i,_(ja,na)), taui, C(_(ic,m),_(jc,n)), work);
    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- unmr3 [real variant] ------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
void
unmr3_impl(Side                      side,
           Transpose                 trans,
           IndexType                 l,
           const GeMatrix<MA>        &A,
           const DenseVector<VTAU>   &tau,
           GeMatrix<MC>              &C,
           DenseVector<VWORK>        &work)
{
    cxxlapack::unmr3(getF77Char(side),
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
                     work.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- unmr3[complex variant] ---------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename MC,
          typename VWORK>
typename RestrictTo<IsInteger<IndexType>::value
                 && IsComplexGeMatrix<MA>::value
                 && IsComplexDenseVector<VTAU>::value
                 && IsComplexGeMatrix<MC>::value
                 && IsComplexDenseVector<VWORK>::value,
         void>::Type
unmr3(Side          side,
      Transpose     trans,
      IndexType     l,
      const MA      &A,
      const VTAU    &tau,
      MC            &&C,
      VWORK         &&work)
{

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType k = tau.length();

    if (side==Left) {
        ASSERT(m>=k);
        ASSERT(A.numCols()==m);
        ASSERT(work.length()==n);
    } else {
        ASSERT(n>=k);
        ASSERT(A.numCols()==n);
        ASSERT(work.length()==m);
    }
    ASSERT(A.numRows()==k);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename RemoveRef<VWORK>::Type VectorWork;

    typename MatrixC::NoView        C_org    = C;
    typename VectorWork::NoView     work_org = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::unmr3_impl(side, trans, l, A, tau, C, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixC::NoView        C_generic    = C;
    typename VectorWork::NoView     work_generic = work;

//
//  restore output arguments
//
    C    = C_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::unmr3_impl(side, trans, l, A, tau, C, work);

    bool failed = false;
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
        std::cerr << "error in: ormr3.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: ormr3.tcc" << std::endl;
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_IMPL_UNMR3_TCC
