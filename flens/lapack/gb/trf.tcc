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

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

/*
namespace generic {

//-- (gb)trf [real and complex variant] ----------------------------------------


template <typename MA, typename VP>
typename GbMatrix<MA>::IndexType
trf_impl(GbMatrix<MA> &A, DenseVector<VP> &piv)
{
    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType  Zero(0), One(1);

    const IndexType    nbMax  = 64;
    const IndexType    ldWork = nbMax+1;
*     ..
*     .. Local Scalars ..
    INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
    DOUBLE PRECISION   TEMP
*     ..
*     .. Local Arrays ..
    DOUBLE PRECISION   WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*   .. External Functions ..
    INTEGER            IDAMAX, ILAENV
    EXTERNAL           IDAMAX, ILAENV
*   ..
*   .. External Subroutines ..
    EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL,
     $                   DSWAP, DTRSM, XERBLA
*   ..
*   .. Intrinsic Functions ..
    INTRINSIC          MAX, MIN
*   ..
*   .. Executable Statements ..
//
//  KV is the number of superdiagonals in the factor U, allowing for
//  fill-in
//
    const IndexType  kl = A.subDiags();
    const IndexType  kv = A.superDiags();
    const IndexType  ku = kv - kl;

    const IndexType  m  = A.numRows();
    const IndexType  n  = A.numCols();
    const IndexType  mn = min(m,n);

    IndexType info = 0;
//
//  Quick return if possible
//
    if (m==0 || n==0) {
        return info;
    }
//
//  Determine the block size for this environment
//
    IndexType nb = ilaenv<T>(1, "GBTRF", "", m, n, kl, ku);
//
//  The block size must not exceed the limit set by the size of the
//  local arrays WORK13 and WORK31.
//
    nb = min(nb, nbMax);

    if (nb<=1 || nb>kl) {
//
//      Use unblocked code
//
        info = tf2(A, piv);
    } else {
//
//      Use blocked code
//
//      Zero the superdiagonal elements of the work array WORK13
//
        Work13.strictUpper() = Zero;
//
//      Zero the subdiagonal elements of the work array WORK31
//
        Work31.strictLower() = Zero;
//
//      Gaussian elimination with partial pivoting
//
//      Set fill-in elements in diagonals KU+2 to KV to zero
//
        for (IndexType d=ku+2; d<=min(kv,n); ++d) {
            A.diag(d) = Zero;
        }
//
//      JU is the index of the last column affected by the current
//      stage of the factorization
//
        IndexType ju = 1;

        for (IndexType j=1; j<=mn; j+=nb) {
            const IndexType jb = min(nb, mn-j+1);
//
//          The active part of the matrix is partitioned
//
//             A11   A12   A13
//             A21   A22   A23
//             A31   A32   A33
//
//          Here A11, A21 and A31 denote the current block of JB columns
//          which is about to be factorized. The number of rows in the
//          partitioning are JB, I2, I3 respectively, and the numbers
//          of columns are JB, J2, J3. The superdiagonal elements of A13
//          and the subdiagonal elements of A31 lie outside the band.
//
            const IndexType i2 = min(kl-jb, m-j-jb+1);
            const IndexType i3 = min(jb, m-j-kl+1);
//
//          J2 and J3 are computed after JU has been updated.
//
//          Factorize the current block of JB columns
//
            DO 80 JJ = J, J + JB - 1
//
//              Set fill-in elements in column JJ+KV to zero
//
                IF( JJ+KV.LE.N ) THEN
                    DO 70 I = 1, KL
                        AB( I, JJ+KV ) = ZERO
   70               CONTINUE
                END IF
//
//              Find pivot and test for singularity. KM is the number of
//              subdiagonal elements in the current column.
//
                KM = MIN( KL, M-JJ )
                JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
                IPIV( JJ ) = JP + JJ - J
                IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                    JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                    IF( JP.NE.1 ) THEN
//
//                      Apply interchange to columns J to J+JB-1
//
                        IF( JP+JJ-1.LT.J+KL ) THEN

                            CALL DSWAP( JB, AB( KV+1+JJ-J, J ),  LDAB-1,
     $                                      AB( KV+JP+JJ-J, J ), LDAB-1 )
                        ELSE
//
//                          The interchange affects columns J to JJ-1 of A31
//                          which are stored in the work array WORK31
//
                            CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                            CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                        END IF
                    END IF
//
//                  Compute multipliers
//
                    CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
//
//                  Update trailing submatrix within the band and within
//                  the current block. JM is the index of the last column
//                  which needs to be updated.
//
                    JM = MIN( JU, J+JB-1 )
                    IF( JM.GT.JJ )
     $                  CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
                ELSE
//
//                  If pivot is zero, set INFO to the index of the pivot
//                  unless a zero pivot has already been found.
//
                    IF( INFO.EQ.0 )
     $                  INFO = JJ
                END IF
//
//              Copy current column of A31 into the work array WORK31
//
                NW = MIN( JJ-J+1, I3 )
                IF( NW.GT.0 )
     $              CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
//
//              Apply the row interchanges to the other blocks.
//
                J2 = MIN( JU-J+1, KV ) - JB
                J3 = MAX( 0, JU-J-KV+1 )
//
//              Use DLASWP to apply the row interchanges to A12, A22, and
//              A32.
//
                CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
//
//              Adjust the pivot indices.
//
                DO 90 I = J, J + JB - 1
                    IPIV( I ) = IPIV( I ) + J - 1
   90           CONTINUE
//
//              Apply the row interchanges to A13, A23, and A33
//              columnwise.
//
                K2 = J - 1 + JB + J2
                DO 110 I = 1, J3
                    JJ = K2 + I
                    DO 100 II = J + I - 1, J + JB - 1
                        IP = IPIV( II )
                        IF( IP.NE.II ) THEN
                            TEMP = AB( KV+1+II-JJ, JJ )
                            AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                            AB( KV+1+IP-JJ, JJ ) = TEMP
                        END IF
  100               CONTINUE
  110           CONTINUE
//
//              Update the relevant part of the trailing submatrix
//
                IF( J2.GT.0 ) THEN
//
//                  Update A12
//
                    CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )

                    IF( I2.GT.0 ) THEN
//
//                      Update A22
//
                        CALL DGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                    END IF

                    IF( I3.GT.0 ) THEN
//
//                      Update A32
//
                        CALL DGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                    END IF
                END IF

                IF( J3.GT.0 ) THEN
//
//                  Copy the lower triangle of A13 into the work array
//                  WORK13
//
                    DO 130 JJ = 1, J3
                        DO 120 II = JJ, JB
                            WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                   CONTINUE
  130               CONTINUE
//
//                  Update A13 in the work array
//
                    CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )

                    IF( I2.GT.0 ) THEN
//
//                      Update A23
//
                        CALL DGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                    END IF

                    IF( I3.GT.0 ) THEN
//
//                      Update A33
//
                        CALL DGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                    END IF
//
//                  Copy the lower triangle of A13 back into place
//
                    DO 150 JJ = 1, J3
                        DO 140 II = JJ, JB
                            AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                   CONTINUE
  150               CONTINUE
                END IF
            ELSE
//
//              Adjust the pivot indices.
//
                DO 160 I = J, J + JB - 1
                    IPIV( I ) = IPIV( I ) + J - 1
  160           CONTINUE
            END IF
//
//          Partially undo the interchanges in the current block to
//          restore the upper triangular form of A31 and copy the upper
//          triangle of A31 back into place
//
            DO 170 JJ = J + JB - 1, J, -1
                JP = IPIV( JJ ) - JJ + 1
                IF( JP.NE.1 ) THEN
//
//                  Apply interchange to columns J to JJ-1
//
                    IF( JP+JJ-1.LT.J+KL ) THEN
//
//                      The interchange does not affect A31
//
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                    ELSE
//
//                      The interchange does affect A31
//
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                    END IF
                END IF
//
//              Copy the current column of A31 back into place
//
                NW = MIN( I3, JJ-J+1 )
                IF( NW.GT.0 )
     $              CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
        }
    }
    END IF

    RETURN
}

} // namespace generic

*/

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

#ifdef USE_CXXLAPACK

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

//
//  Call implementation
//
    //IndexType info = LAPACK_SELECT::trf_impl(A, piv);
    IndexType info = external::trf_impl(A, piv);

    return info;
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GB_TRF_TCC
