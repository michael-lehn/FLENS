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
       SUBROUTINE DLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LARZ_TCC
#define FLENS_LAPACK_LA_LARZ_TCC 1

#include <flens/lapack/lapack.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- larz [real variant] -------------------------------------------------------

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larz_impl(Side                   side,
          const DenseVector<VV>  &v,
          const TAU              &tau,
          GeMatrix<MC>           &C,
          DenseVector<VWORK>     &work)
{
    typedef typename GeMatrix<MC>::ElementType  ElementType;
    typedef typename GeMatrix<MC>::IndexType    IndexType;

    const Underscore<IndexType>  _;

    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType l = v.length();

    const ElementType  Zero(0), One(1);

    if (side==Left) {
//
//      Form  H * C
//
        if (tau!=Zero) {
//
//          w( 1:n ) = C( 1, 1:n )
//
            work = C(1,_);
//
//          w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )
//
            blas::mv(Trans, One, C(_(m-l+1,m),_), v, One, work);
//
//          C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
//
            C(1,_) -= tau*work;
//
//          C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
//                              tau * v( 1:l ) * w( 1:n )**T
//
            blas::r(-tau, v, work, C(_(m-l+1,m),_));
        }

    } else {
//
//      Form  C * H
//
        if (tau!=Zero) {
//
//          w( 1:m ) = C( 1:m, 1 )
//
            work = C(_,1);
//
//          w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
//
            blas::mv(NoTrans, One, C(_,_(n-l+1,n)), v, One, work);
//
//          C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
//
            C(_,1) -= tau*work;
//
//          C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
//                              tau * w( 1:m ) * v( 1:l )**T
//
            blas::r(-tau, work, v, C(_,_(n-l+1,n)));

        }

    }
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- larz [real variant] -------------------------------------------------------

template <typename VV, typename TAU, typename MC, typename VWORK>
void
larz_impl(Side                   side,
          const DenseVector<VV>  &v,
          const TAU              &tau,
          GeMatrix<MC>           &C,
          DenseVector<VWORK>     &work)
{
    typedef typename GeMatrix<MC>::IndexType  IndexType;

    cxxlapack::larz<IndexType>(getF77Char(side),
                               C.numRows(),
                               C.numCols(),
                               v.length(),
                               v.data(),
                               v.inc()*v.stride(),
                               tau,
                               C.data(),
                               C.leadingDimension(),
                               work.data());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- larz [real variant] -------------------------------------------------------

template <typename VV, typename TAU, typename MC, typename VWORK>
typename RestrictTo<IsRealDenseVector<VV>::value
                 && IsReal<TAU>::value
                 && IsRealGeMatrix<MC>::value
                 && IsRealDenseVector<VWORK>::value,
         void>::Type
larz(Side       side,
     const VV   &v,
     const TAU  &tau,
     MC         &&C,
     VWORK      &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MC>::Type    MatrixC;
    typedef typename RemoveRef<VWORK>::Type VectorWork;
    typedef typename MatrixC::IndexType     IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG
    const IndexType m = C.numRows();
    const IndexType n = C.numCols();
    const IndexType l = v.length();

    if (side==Left) {
        ASSERT(m>=l);
        ASSERT(work.length()==n);
    } else {
        ASSERT(n>=l);
        ASSERT(work.length()==m);
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixC::NoView        C_org    = C;
    typename VectorWork::NoView     work_org = work;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::larz_impl(side, v, tau, C, work);

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
    external::larz_impl(side, v, tau, C, work);

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
        std::cerr << "error in: larz.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: larz.tcc" << std::endl;
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LARZ_TCC
