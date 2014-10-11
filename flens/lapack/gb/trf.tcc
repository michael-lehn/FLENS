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
       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
       SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GB_TRF_TCC
#define FLENS_LAPACK_GB_TRF_TCC 1

#include <cxxstd/algorithm.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

//-- (gb)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GbMatrix<MA>::IndexType
trf_impl(GbMatrix<MA> &A, DenseVector<VP> &piv)
{
    using std::min;
    using std::max;

    typedef typename GbMatrix<MA>::GeNoView        GeNoView;
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
    const IndexType  ldAB = A.leadingDimension();

    const ElementType zero(0), one(1);


    const Underscore<IndexType> _;

    IndexType nbmax = 64, ldWork = nbmax+1;



    IndexType nb = lapack::ilaenv<ElementType>(1, "GBTRF", "", m, n, kl, ku);

//
//  The block size must not exceed the limit set by the size of the
//  local arrays Work13 and Work31.
//
    nb = min( nb, nbmax );

    if( nb<=1 || nb>kl ) {
//
//      Use unblocked code
//
        return lapack::tf2(A, piv);

    } else {

        GeNoView Work13(ldWork, nbmax), Work31(ldWork, nbmax);
//
//      Use blocked code
//

//
//      Gaussian elimination with partial pivoting
//
//      Set fill-in elements in columns ku+2 to kv to zero
//
        for(IndexType j = ku + 2; j<=min( kv, n ); ++j) {
            AB(_(kv - j + 2, kl), j) = zero;
        }
//
//      ju is the index of the last column affected by the current
//      stage of the factorization
//
        IndexType ju = 1;

        for (IndexType j = 1; j<=min( m, n ); j+=nb) {
            IndexType jb = min( nb, min( m, n )-j+1 );
//
//          The active part of the matrix is partitioned
//
//              A11   A12   A13
//              A21   A22   A23
//              A31   A32   A33
//
//          Here A11, A21 and A31 denote the current block of jB columns
//          which is about to be factorized. The number of rows in the
//          partitioning are jb, i2, i3 respectively, and the numbers
//          of columns are jB, j2, Jj3. The superdiagonal elements of A13
//          and the subdiagonal elements of A31 lie outside the band.
//
            IndexType i2 = min( kl-jb, m-j-jb+1 );
            IndexType i3 = min( jb, m-j-kl+1 );
//
//          j2 and j3 are computed after ju has been updated.
//
//          Factorize the current block of jb columns
//
            for (IndexType jj = j; jj<=j + jb - 1; ++jj) {
//
//              Set fill-in elements in column JJ+KV to zero
//
                if( jj+kv<=n ) {
                    AB(_(1,kl),jj+kv) = zero;
                }
//
//              Find pivot and test for singularity. km is the number of
//              subdiagonal elements in the current column.
//
                IndexType km = min( kl, m-jj );
                IndexType jp = blas::iamax(AB(_(kv+1,kv+km+1),jj));
                piv( jj ) = jp + jj - j;
                if( AB( kv+jp, jj )!=zero ) {
                    ju = max( ju, min( jj+ku+jp-1, n ) );
                    if( jp!=1 ) {
//
//                      Apply interchange to columns j to j+jb-1
//
                        if( jp+jj-1<j+kl ) {

                            blas::swap(A(jj,_(j,j+jb-1)),
                                       A(jj+jp-1,_(j,j+jb-1)));

                        } else {
//
//                          The interchange affects columns j to jj-1 of A31
//                          which are stored in the work array WORK31
//
                            blas::swap(A(jj,_(j,jj-1)),
                                       Work31(jp+jj-j-kl, _(1, jj-j)));
                            blas::swap(A(jj,_(jj,j+jb-1)),
                                       A(jj+jp-1,_(jj,j+jb-1)));

                        }
                    }
//
//                  Compute multipliers
//
                    blas::scal(one / AB( kv+1, jj ), AB( _(kv+2,kv+km+1), jj ));

//
//                  Update trailing submatrix within the band and within
//                  the current block. jm is the index of the last column
//                  which needs to be updated.
//
                    IndexType jm = min( ju, j+jb-1 );
                    if( jm>jj ) {
                        GeMatrixView AB_tmp = FullStorageView(km, jm-jj,
                                                              &AB( kv+1, jj+1),
                                                              ldAB-1);

                        blas::r(-one, AB( _(kv+2, kv+km+1), jj ),
                                A(jj,_(jj+1,jm)), AB_tmp);

                    }
                } else {
//
//                  If pivot is zero, return index
//
                    return jj;
                }
//
//              Copy current column of A31 into the work array WORK31
//
                IndexType NW = min( jj-j+1, i3 );
                if( NW>0 ) {
                    blas::copy(AB( _(kv+kl+1-jj+j,kv+kl-jj+j+NW), jj ),
                               Work31( _(1,NW), jj-j+1 ));
                }
            }
            if( j+jb<=n ) {
//
//              Apply the row interchanges to the other blocks.
//
                IndexType j2 = min( ju-j+1, kv ) - jb;
                IndexType j3 = max( 0, ju-j-kv+1 );
//
//              Use DLASWP to apply the row interchanges to A12, A22, and
//              A32.
//
                const IndexType mAB_tmp = min(ldAB-1, AB.numRows()-kv+jb);
                GeMatrixView AB_tmp = FullStorageView(mAB_tmp, j2,
                                                      &AB( kv+1-jb, j+jb ),
                                                      ldAB-1);
                lapack::laswp(AB_tmp, piv(_(j,j+jb-1)));

//
//              Adjust the pivot indices.
//
                for(IndexType i = j; i<=j + jb - 1; ++i) {
                    piv( i ) = piv( i ) + j - 1;
                }
//
//              Apply the row interchanges to A13, A23, and A33
//              columnwise.
//
                IndexType k2 = j - 1 + jb + j2;
                for (IndexType i = 1; i<=j3; ++i) {
                    IndexType jj = k2 + i;
                    for (IndexType ii = j + i - 1; ii<=j + jb - 1; ++ii) {
                        IndexType ip = piv( ii );
                        if ( ip!=ii ) {
                            ElementType TEMP     = AB( kv+1+ii-jj, jj );
                            AB( kv+1+ii-jj, jj ) = AB( kv+1+ip-jj, jj );
                            AB( kv+1+ip-jj, jj ) = TEMP;
                        }
                    }
                }
//
//              Update the relevant part of the trailing submatrix
//
                if( j2>0 ) {
//
//                  Update A12
//
                    GeMatrixView AB_tmp_1 = FullStorageView(jb, jb,
                                                            &AB(kv+1, j),
                                                            ldAB-1);
                    GeMatrixView AB_tmp_2 = FullStorageView(jb, j2,
                                                            &AB(kv+1-jb, j+jb),
                                                            ldAB-1);

                    blas::sm(Left, NoTrans, one, AB_tmp_1.lowerUnit(),
                             AB_tmp_2);

                    if( i2>0 ) {
//
//                      Update A22
//
                        GeMatrixView AB_tmp_1 = FullStorageView(i2, jb,
                                                        &AB(kv+1+jb, j),
                                                        ldAB-1);
                        GeMatrixView AB_tmp_2 = FullStorageView(jb, j2,
                                                        &AB( kv+1-jb, j+jb),
                                                        ldAB-1);
                        GeMatrixView AB_tmp_3 = FullStorageView(i2, j2,
                                                        &AB( kv+1, j+jb),
                                                        ldAB-1);
                        AB_tmp_3 -= AB_tmp_1*AB_tmp_2;


                    }

                    if( i3>0 ) {
//
//                      Update A32
//
                        GeMatrixView AB_tmp_1 = FullStorageView(jb, j2,
                                                        &AB(kv+1-jb, j+jb),
                                                        ldAB-1);
                        GeMatrixView AB_tmp_2 = FullStorageView(i3, j2,
                                                        &AB( kv+kl+1-jb, j+jb),
                                                        ldAB-1);
                        AB_tmp_2 -= Work31(_(1,i3),_(1,jb))*AB_tmp_1;

                    }

                }

                if( j3>0 ) {
//
//                  Copy the lower triangle of A13 into the work array
//                  WORK13
//
                    for (IndexType jj=1; jj<=j3; ++jj) {
                        Work13(_(jj,jb), jj) = AB(_(1,jb-jj+1), jj+j+kv-1);
                    }

//
//                  Update A13 in the work array
//
                    GeMatrixView AB_tmp = FullStorageView(jb, jb, &AB(kv+1, j),
                                                          ldAB-1);


                    blas::sm(Left, NoTrans, one, AB_tmp.lowerUnit(),
                             Work13(_(1,jb),_(1,j3)));

                    if( i2>0 ) {
//
//                      Update A23
//
                        GeMatrixView AB_tmp_1 = FullStorageView(i2, jb,
                                                                &AB(kv+1+jb, j),
                                                                ldAB-1);
                        GeMatrixView AB_tmp_3 = FullStorageView(i2, j3,
                                                                &AB(1+jb, j+kv),
                                                                ldAB-1);
                        AB_tmp_3 -= AB_tmp_1*Work13(_(1,jb),_(1,j3));

                    }

                    if( i3>0 ) {
//
//                      Update A33
//
                        GeMatrixView AB_tmp = FullStorageView(i3, j3,
                                                              &AB( 1+kl, j+kv),
                                                              ldAB-1);
                        AB_tmp -= Work31(_(1,i3),_(1,jb))
                                 *Work13(_(1,jb),_(1,j3));
                    }
//
//                  Copy the lower triangle of A13 back into place
//
                    for (IndexType jj = 1; jj<=j3; ++jj) {
                        AB(_(1,jb-jj+1), jj+j+kv-1) = Work13(_(jj,jb),jj) ;
                    }
                }
            } else {
//
//              Adjust the pivot indices.
//
                for (IndexType i = j; i<= j + jb - 1; ++i) {
                    piv( i ) = piv( i ) + j - 1;
                }
            }
//
//          Partially undo the interchanges in the current block to
//          restore the upper triangular form of A31 and copy the upper
//          triangle of A31 back into place
//
            for (IndexType jj = j + jb - 1; jj>= j; --jj) {
                IndexType jp = piv( jj ) - jj + 1;
                if( jp!=1 ) {
//
//                  Apply interchange to columns j to jj-1
//
                    if( jp+jj-1<j+kl ) {
//
//                      The interchange does not affect A31
//
                        blas::swap(A(jj,_(j,jj-1)),A(jj+jp-1,_(j,jj-1)));

                    } else {
//
//                      The interchange does affect A31
//
                        blas::swap(A(jj,_(j,jj-1)), Work31(jp+jj-j-kl, _(1,jj-j)));

                    }
                }
//
//              Copy the current column of A31 back into place
//
                IndexType nw = min( i3, jj-j+1 );
                if( nw>0 ) {
                    blas::copy(Work31(_(1,nw), jj-j+1),
                               AB( _(kv+kl+1-jj+j,kv+kl-jj+j+nw), jj ));
                }
            }
        }
    }
    return 0;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (gb)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VP>
typename GbMatrix<MA>::IndexType
trf_impl(GbMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    const IndexType  k = A.numSuperDiags()-A.numSubDiags();
    return cxxlapack::gbtrf<IndexType>(A.numRows(), A.numCols(),
                                       A.numSubDiags(), k,
                                       A.data(), A.leadingDimension(),
                                       piv.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

//-- (gb)trf [real and complex variant] ----------------------------------------

template <typename MA, typename VPIV>
typename RestrictTo<IsGbMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
trf(MA &&A, VPIV &&piv)
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
    IndexType info = LAPACK_SELECT::trf_impl(A, piv);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    IndexType _info = external::trf_impl(_A, _piv);

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

#endif // FLENS_LAPACK_GB_TRF_TCC
