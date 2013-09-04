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
       SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LATRZ_TCC
#define FLENS_LAPACK_LA_LATRZ_TCC 1

#include <flens/lapack/lapack.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- latrz [real variant] ------------------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename
RestrictTo<IsRealGeMatrix<GeMatrix<MA> >::value,
           void>::Type
latrz_impl(IndexType             l,
           GeMatrix<MA>          &A,
           DenseVector<VTAU>     &tau,
           DenseVector<VWORK>    &work)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;

    const Underscore<IndexType>  _;

    const ElementType  Zero(0);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//  Quick return if possible
//
    if (m==0) {
        return;
    } else if (m==n) {
        tau = Zero;
        return;
    }

    for (IndexType i=m; i>=1; --i) {
//
//      Generate elementary reflector H(i) to annihilate
//      [ A(i,i) A(i,n-l+1:n) ]
//
        larfg(l+1, A(i,i), A(i,_(n-l+1,n)), tau(i));

//
//      Apply H(i) to A(1:i-1,i:n) from the right
//
        auto _work = work(_(1,i-1));
        larz(Right, A(i,_(n-l+1,n)), tau(i), A(_(1,i-1),_(i,n)), _work);
    }
}

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename
RestrictTo<IsComplexGeMatrix<GeMatrix<MA> >::value,
           void>::Type
latrz_impl(IndexType             l,
           GeMatrix<MA>          &A,
           DenseVector<VTAU>     &tau,
           DenseVector<VWORK>    &work)
{
    
    typedef typename GeMatrix<MA>::ElementType  ElementType;

    const Underscore<IndexType>  _;

    const ElementType  Zero(0);
    ElementType alpha(0);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

//
//  Quick return if possible
//
    if (m==0) {
        return;
    } else if (m==n) {
        tau = Zero;
        return;
    }

    for (IndexType i=m; i>=1; --i) {
//
//      Generate elementary reflector H(i) to annihilate
//      [ A(i,i) A(i,n-l+1:n) ]
//

        A(i, _(n-l+1,n)) = conjugate(A(i, _(n-l+1,n)));
        alpha = cxxblas::conjugate(A(i,i));
        larfg(l+1, alpha, A(i,_(n-l+1,n)), tau(i));
        tau(i) = cxxblas::conjugate(tau(i));
//
//      Apply H(i) to A(1:i-1,i:n) from the right
//
        auto _work = work(_(1,i-1));
        larz(Right, A(i,_(n-l+1,n)), cxxblas::conjugate(tau(i)), A(_(1,i-1),_(i,n)), _work);
        A(i,i) = cxxblas::conjugate(alpha);
    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- latrz [real and complex variant] ------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
void
latrz_impl(IndexType             l,
           GeMatrix<MA>          &A,
           DenseVector<VTAU>     &tau,
           DenseVector<VWORK>    &work)
{
    cxxlapack::latrz(A.numRows(),
                     A.numCols(),
                     l,
                     A.data(),
                     A.leadingDimension(),
                     tau.data(),
                     work.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- latrz[real and complex variant] -------------------------------------------

template <typename IndexType, typename MA, typename VTAU, typename VWORK>
typename RestrictTo<IsInteger<IndexType>::value
                 && IsGeMatrix<MA>::value
                 && IsDenseVector<VTAU>::value
                 && IsDenseVector<VWORK>::value,
         void>::Type
latrz(IndexType     l,
      MA            &&A,
      VTAU          &&tau,
      VWORK         &&work)
{

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = A.numRows();
    const IndexType n = A.numCols();
    ASSERT(n-m>=l);
    ASSERT(tau.length()==m);
    ASSERT(work.length()==m);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VTAU>::Type  VectorTau;
    typedef typename RemoveRef<VWORK>::Type VectorWork;

    typename MatrixA::NoView        A_org    = A;
    typename VectorTau::NoView      tau_org  = tau;
    typename VectorWork::NoView     work_org = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::latrz_impl(l, A, tau, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixA::NoView        A_generic    = A;
    typename VectorTau::NoView      tau_generic  = tau;
    typename VectorWork::NoView     work_generic = work;

//
//  restore output arguments
//
    A    = A_org;
    tau  = tau_org;
    work = work_org;

//
//  Compare generic results with results from the native implementation
//
    external::latrz_impl(l, A, tau, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(tau_generic, tau, "tau_generic", "tau")) {
        std::cerr << "CXXLAPACK: tau_generic = " << tau_generic << std::endl;
        std::cerr << "F77LAPACK: tau = " << tau << std::endl;
        failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: latrz.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: latrz.tcc" << std::endl;
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LATRZ_TCC
