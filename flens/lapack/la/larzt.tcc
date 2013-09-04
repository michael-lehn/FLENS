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
       SUBROUTINE DLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
       SUBROUTINE ZLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_LA_LARZT_TCC
#define FLENS_LAPACK_LA_LARZT_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- larzt [real variant] ------------------------------------------------------

template <typename MV, typename VTAU, typename MT>
typename 
RestrictTo<IsRealGeMatrix<GeMatrix<MV> >::value,
           void>::Type
larzt_impl(Direction                 direction,
           StoreVectors              storeVectors,
           GeMatrix<MV>              &V,
           const DenseVector<VTAU>   &tau,
           TrMatrix<MT>              &T)
{
    typedef typename GeMatrix<MV>::ElementType  ElementType;
    typedef typename GeMatrix<MV>::IndexType    IndexType;

    const Underscore<IndexType>  _;

    const ElementType  Zero(0);

//
//  Check for currently supported options
//
    ASSERT(storeVectors==RowWise);
    ASSERT(direction==Backward);

    const IndexType k = T.dim();

    for (IndexType i=k; i>=1; --i) {
        if (tau(i)==Zero) {
//
//          H(i)  =  I
//
            T(_(i,k),i) = Zero;
        } else {
//
//          general case
//
            if (i<k) {
//
//              T(i+1:k,i) = - tau(i) * V(i+1:k,1:n)**T * V(i,1:n)
//
                const auto rows = _(i+1,k);
                blas::mv(NoTrans, -tau(i), V(rows,_), V(i,_), Zero, T(rows,i));
//
//              T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
//
                blas::mv(NoTrans, T(rows,rows).lower(), T(rows,i));
            }
            T(i,i) = tau(i);
        }
    }
}

//-- larzt [complex variant] ---------------------------------------------------

template <typename MV, typename VTAU, typename MT>
typename 
RestrictTo<IsComplexGeMatrix<GeMatrix<MV> >::value,
           void>::Type
larzt_impl(Direction                 direction,
           StoreVectors              storeVectors,
           GeMatrix<MV>              &V,
           const DenseVector<VTAU>   &tau,
           TrMatrix<MT>              &T)
{
    typedef typename GeMatrix<MV>::ElementType  ElementType;
    typedef typename GeMatrix<MV>::IndexType    IndexType;

    const Underscore<IndexType>  _;

    const ElementType  Zero(0);

//
//  Check for currently supported options
//
    ASSERT(storeVectors==RowWise);
    ASSERT(direction==Backward);

    const IndexType k = T.dim();
    const IndexType n = V.numCols();
    
    for (IndexType i=k; i>=1; --i) {
        if (tau(i)==Zero) {
//
//          H(i)  =  I
//
            T(_(i,k),i) = Zero;
        } else {
//
//          general case
//
            if (i<k) {
//
//              T(i+1:k,i) = - tau(i) * V(i+1:k,1:n)**H * V(i,1:n)
//
                const auto rows = _(i+1,k);
                V(i, _(1, n)) = conjugate(V(i, _(1, n)));
                blas::mv(NoTrans, -tau(i), V(rows,_), V(i,_), Zero, T(rows,i));
                V(i, _(1, n)) = conjugate(V(i, _(1, n)));
//
//              T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
//
                blas::mv(NoTrans, T(rows,rows).lower(), T(rows,i));
            }
            T(i,i) = tau(i);
        }
    }
}


} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- larzt [real variant] ------------------------------------------------------

template <typename MV, typename VTAU, typename MT>
void
larzt_impl(Direction                 direction,
           StoreVectors              storeVectors,
           GeMatrix<MV>              &V,
           const DenseVector<VTAU>   &tau,
           TrMatrix<MT>              &T)
{
    typedef typename GeMatrix<MV>::IndexType  IndexType;
    const IndexType n = (storeVectors==RowWise) ? V.numCols()
                                                : V.numRows();
    cxxlapack::larzt(getF77Char(direction),
                     getF77Char(storeVectors),
                     n,
                     T.dim(),
                     V.data(),
                     V.leadingDimension(),
                     tau.data(),
                     T.data(),
                     T.leadingDimension());
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- larzt[real and complex variant] -------------------------------------------

template <typename MV, typename VTAU, typename MT>
typename RestrictTo<IsGeMatrix<MV>::value
                 && IsDenseVector<VTAU>::value
                 && IsTrMatrix<MT>::value,
         void>::Type
larzt(Direction      direction,
      StoreVectors   storeVectors,
      MV             &&V,
      const VTAU     &tau,
      MT             &&T)
{

//
//  Test the input parameters
//
#   ifndef NDEBUG

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MV>::Type    MatrixV;
    typedef typename MatrixV::IndexType     IndexType;


//  From LAPACK: Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
    ASSERT(storeVectors==RowWise);
    ASSERT(direction==Backward);

    const IndexType k = T.dim();
    const IndexType n = (storeVectors==RowWise) ? V.numCols()
                                                : V.numRows();
    ASSERT(tau.length()==k);
    if (storeVectors==RowWise) {
        ASSERT(V.numRows()==k);
    } else {
        ASSERT(V.numRows()==n);
    }
    ASSERT(T.upLo()==Lower);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<MT>::Type    MatrixT;
    
    typename MatrixV::NoView        V_org = V;
    typename MatrixT::NoView        T_org = T;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::larzt_impl(direction, storeVectors, V, tau, T);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename MatrixV::NoView        V_generic = V;
    typename MatrixT::NoView        T_generic = T;

//
//  restore output arguments
//
    V = V_org;
    T = T_org;

//
//  Compare generic results with results from the native implementation
//
    external::larzt_impl(direction, storeVectors, V, tau, T);

    bool failed = false;
    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }
    if (! isIdentical(T_generic, T, "T_generic", "T")) {
        std::cerr << "CXXLAPACK: T_generic = " << T_generic << std::endl;
        std::cerr << "F77LAPACK: T = " << T << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: larzt.tcc" << std::endl;
        ASSERT(0);
    } else {
        // std::cerr << "passed: larzt.tcc" << std::endl;
    }
#   endif

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LARZT_TCC
