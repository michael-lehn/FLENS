/*
 *   Copyright (c) 2013, Michael Lehn, Klaus Pototzky
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
       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
       SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GB_TF2_TCC
#define FLENS_LAPACK_GB_TF2_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (gb)tf2 [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GbMatrix<MA>::IndexType
tf2_impl(GbMatrix<MA> &A, DenseVector<VP> &piv)
{
	using std::min;
	using std::max;

	typedef typename GbMatrix<MA>::GeView          GeMatrixView;
	typedef typename GbMatrix<MA>::FullStorageView FullStorageView;
    typedef typename GbMatrix<MA>::IndexType       IndexType;
    typedef typename GbMatrix<MA>::ElementType     ElementType;

    GeMatrixView AB = A.viewStorageGeMatrix();

    const IndexType  m  = A.numRows();
    const IndexType  n  = A.numCols();
    const IndexType  kl = A.numSubDiags();
    const IndexType  ku = A.numSuperDiags()-A.numSubDiags();
    const IndexType  kv = kl + ku;
    const IndexType  ldAB = AB.leadingDimension();

    const ElementType zero(0), one(1);


    const Underscore<IndexType> _;

//
//  Gaussian elimination with partial pivoting
//
//  Set fill-in elements in columns ku+2 to kv to zero.
//
	for( IndexType j = ku + 2; j<=min( kv, n ); ++j) {
		AB(_( kv - j + 2, kl), j) = zero;
	}
//
//  ju is the index of the last column affected by the current stage
//  of the factorization.
//
	IndexType ju = 1;

	for (IndexType j = 1; j<= min( m, n ); ++j) {
//
//      Set fill-in elements in column j+kv to zero.
//
		if( j+kv<=n ) {;
			AB(_(1,kl),j+kv) = zero;
		}
//
//        Find pivot and test for singularity. km is the number of
//        subdiagonal elements in the current column.
//
		IndexType km = min( kl, m-j );
		IndexType jp = blas::iamax(AB( _(kv+1,kv+km+1), j ));
		piv( j ) = jp + j - 1;
		if( AB( kv+jp, j ) != zero ) {
			ju = max( ju, min( j+ku+jp-1, n ) );
//
//          Apply interchange to columns j to ju.
//
			if( jp!=1 ) {
				blas::swap(A(jp+j-1,_(j,ju)), A(j,_(j,ju)));
			}

			if( km>0 ) {
//
//              Compute multipliers.
//
				blas::scal(one / AB( kv+1, j ), AB(_(kv+2, kv+km+1), j ));
//
//              Update trailing submatrix within the band.
//
				if( ju>j ) {
					GeMatrixView  AB_tmp = FullStorageView( km, ju-j, &AB( kv+1, j+1 ), ldAB-1);
					blas::r(-one, AB(_(kv+2,kv+km+1),j), A(j,_(j+1,ju)), AB_tmp);

				}
			}
		} else {
//
//          If pivot is zero, return the index
//
			return j;
		}
	}
	return 0;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (gb)tf2 [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GbMatrix<MA>::IndexType
tf2_impl(GbMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    const IndexType  k = A.numSuperDiags()-A.numSubDiags();
    return cxxlapack::gbtf2<IndexType>(A.numRows(), A.numCols(),
                                       A.numSubDiags(), k,
                                       A.data(), A.leadingDimension(),
                                       piv.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (gb)tf2 [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
tf2(MA &&A, VPIV &&piv)
{
    using std::min;

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    const IndexType  mn = min(A.numRows(), A.numCols());

    if (piv.length()!=mn) {
        piv.resize(mn);
    }
    ASSERT(piv.length()==mn);

#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT((piv.inc()>0 && piv.firstIndex()==1)
        || (piv.inc()<0 && piv.firstIndex()==A.numRows()));
#   endif

#   ifdef CHECK_CXXLAPACK

    typedef typename RemoveRef<VPIV>::Type  VectorPiv;
    
//
//  Make copies of output arguments
//
    typename MatrixA::NoView    _A      = A;
    typename VectorPiv::NoView  _piv    = piv;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::tf2_impl(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::tf2_impl(_A, _piv);

    bool failed = false;
    if (! isIdentical(A, _A, " A", "_A")) {
        std::cerr << "CXXLAPACK:  A = " << A << std::endl;
        std::cerr << "F77LAPACK: _A = " << _A << std::endl;
        failed = true;
    }

    if (! isIdentical(piv, _piv, " piv", "_piv")) {
        std::cerr << "CXXLAPACK:  piv = " << piv << std::endl;
        std::cerr << "F77LAPACK: _piv = " << _piv << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }

#   endif

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GB_TF2_TCC
