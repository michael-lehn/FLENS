/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
      $                   INFO )
       SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
      $                   INFO )
 *
 *  -- LAPACK routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 *
 */

#ifndef FLENS_LAPACK_GB_TRS_TCC
#define FLENS_LAPACK_GB_TRS_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (gb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const GbMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
	using std::min;
	using std::max;

	typedef typename GbMatrix<MA>::ConstGeView     ConstGeMatrixView;
    typedef typename GbMatrix<MA>::IndexType       IndexType;
    typedef typename GbMatrix<MA>::ElementType     ElementType;

    const ConstGeMatrixView AB = A.viewStorageGeMatrix();

    const IndexType  n  = A.numCols();
    const IndexType  kl = A.numSubDiags();
    const IndexType  ku = A.numSuperDiags()-A.numSubDiags();
    const IndexType  kd = kl + ku + 1;
    const IndexType  nrhs = B.numCols();

    ElementType one(1);
    const Underscore<IndexType> _;

    if (trans==NoTrans) {
//
//      Solve  A*X = B.
//
//      Solve L*X = B, overwriting B with X.
//
//      L is represented as a product of permutations and unit lower
//      triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
//      where each transformation L(i) is a rank-one modification of
//      the identity matrix.
//
		if( kl>0 ) {
		   for (IndexType j = 1; j<= n - 1; ++j) {
			   IndexType lm = min( kl, n-j );
			   IndexType l = piv( j );
			  if( l!=j ) {
				  blas::swap(B(l,_), B(j,_));
			  }
			  blas::r(-one, AB(_(kd+1,kd+lm),j), B(j,_), B(_(j+1,j+lm),_));

		   }
		}
		for (IndexType I = 1; I<=nrhs; ++I) {
//
//          Solve U*X = B, overwriting B with X.
//
			blas::sv(NoTrans, A.upper(), B(_,I));
		}
    } else if (trans==Trans){
//
//      Solve A**T*X = B.
//
		for (IndexType i = 1; i<=nrhs; ++i) {
			blas::sv(Trans, A.upper(), B(_,i));
		}
//
//      Solve L**T*X = B, overwriting B with X.
//
		if( kl>0 ) {
			for (IndexType j = n - 1; j>=1; --j) {
				IndexType lm = min( kl, n-j );

				B(j,_) -= transpose(B( _(j+1,j+lm), _(1,nrhs) ))*AB( _(kd+1,kd+lm), j );
				IndexType l = piv( j );
				if( l!=j ) {
					blas::swap(B(l,_),B(j,_));
				}
			}
		}
    } else if (trans==ConjTrans) {
//
//      Solve A**H*X = B.
//
		for (IndexType i = 1; i<=nrhs; ++i) {
			blas::sv(ConjTrans, A.upper(), B(_,i));
		}
//
//      Solve L**T*X = B, overwriting B with X.
//
		if( kl>0 ) {
			for (IndexType j = n - 1; j>=1; --j) {
				IndexType lm = min( kl, n-j );
				B(j,_)  = conjugate(B(j,_));
				B(j,_) -= conjTrans(B( _(j+1,j+lm), _(1,nrhs) ))*AB( _(kd+1,kd+lm), j );
				B(j,_)  = conjugate(B(j,_));
				IndexType l = piv( j );
				if( l!=j ) {
					blas::swap(B(l,_),B(j,_));
				}
			}
		}
    } else {
    	ASSERT(0);
    }
}

}

//== interface for external lapack =============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (gb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VP, typename MB>
void
trs_impl(Transpose trans, const GbMatrix<MA> &A, const DenseVector<VP> &piv,
         GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    IndexType info;
    info = cxxlapack::gbtrs<IndexType>(getF77Char(trans),
                                       A.numRows(),
                                       A.numSubDiags(),
                                       A.numSuperDiags()-A.numSubDiags(),
                                       B.numCols(),
                                       A.data(),
                                       A.leadingDimension(),
                                       piv.data(),
                                       B.data(),
                                       B.leadingDimension());
    ASSERT(info==0);
}

} // namespace external

#endif // USE_CXXLAPACK


//== public interface ==========================================================

//-- (gb)trs [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV, typename MB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, MB &&B)
{
    LAPACK_DEBUG_OUT("(gb)trs [real/complex]");

//
//  Test the input parameters
//
#   ifndef NDEBUG

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename MatrixA::IndexType      IndexType;
    
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());

    const IndexType n = A.numRows();

    ASSERT(piv.firstIndex()==1);
    ASSERT(piv.length()==n);

    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==n);
#   endif

#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<MB>::Type     MatrixB;
    
//
//  Make copies of output arguments
//
    typename MatrixB::NoView  B_org   = B;
#   endif
    

//
//  Call implementation
//
    LAPACK_SELECT::trs_impl(trans, A, piv, B);
    
//
//  Compare results
//
#   ifdef CHECK_CXXLAPACK
    typename MatrixB::NoView  B_generic   = B;

    B   = B_org;

    external::trs_impl(trans, A, piv, B);

    bool failed = false;
    if (! isIdentical(B_generic, B, "B_generic", "B")) {
        std::cerr << "CXXLAPACK: B_generic = " << B_generic << std::endl;
        std::cerr << "F77LAPACK: B = " << B << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    } else {
        // std::cerr << "passed: (ge)trs.tcc" << std::endl;
    }
#   endif

}

//-- (gb)trs [variant if rhs is vector] ----------------------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsDenseVector<VB>::value,
         void>::Type
trs(Transpose trans, const MA &A, const VPIV &piv, VB &&b)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type     MatrixA;
    typedef typename RemoveRef<VB>::Type     VectorB;

//
//  Create matrix view from vector b and call above variant
//
    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    trs(trans, A, piv, B);
}


} } // namespace lapack, flens

#endif // FLENS_LAPACK_GB_TRS_TCC
