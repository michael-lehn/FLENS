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
       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GESV_TRI_TCC
#define FLENS_LAPACK_GESV_TRI_TCC 1

#include <algorithm>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

#include <flens/lapack/interface/include/f77lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri_generic(GeMatrix<MA> &A, DenseVector<VP> &piv, DenseVector<VWORK> &work)
{
    using std::max;

    IndexType info = 0;
    IndexType nb = ilaenv<ElementType>(1, "GETRI", "", n);

    const IndexType lWorkOpt = n*nb;

    if (work.length()==0) {
        work.resize(max(lWorkOpt, IndexType(1)));
    }
    work(1) = lWorkOpt;

//
//  Quick return if possible
//
    if (n==0) {
        return info;
    }
//
//  Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
//  and the inverse is not computed.
//
    info = tri(A.upper());
    if (info>0) {
        return info;
    }

    IndexType nbMin = 2;
    const IndexType ldWork = n;

    IndexType iws;

    if (nb>1 && nb<n) {
        iws = max(ldWork*nb, IndexType(1));
        if (work.length()<iws) {
            nb = lWork / ldWork;
            nbMin = max(IndexType(2), ilaenv<ElementType>(2, "GETRI", "", n));
        }
    } else {
        iws = n;
    }
//
//  Solve the equation inv(A)*L = inv(U) for inv(A).
//
    if (nb<nbMin || nb>=n) {
//
//      Use unblocked code.
//
        DO 20 J = N, 1, -1
//
//          Copy current column of L to WORK and replace with zeros.
//
            DO 10 I = J + 1, N
                WORK( I ) = A( I, J )
                A( I, J ) = ZERO
   10       CONTINUE
//
//          Compute current column of inv(A).
//
            IF( J.LT.N )
     $          CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20   CONTINUE
    } else {
//
//      Use blocked code.
//
        const IndexType nn = ((n-1)/nb)*nb + 1;
        for (IndexType j=nn; j>=1; j-=nb) {
            const IndexType jb = min(nb, n-j+1);
//
//          Copy current block column of L to WORK and replace with
//          zeros.
//
            for (IndexType jj=j; jj<=j+nb-1; ++jj) {
                for (IndexType i=jj+1; i<=n; ++i) {
                    Work(i,jj-j) = A(i,jj);
                    A(i, jj) = Zero;
                }
            }
//
//          Compute current block column of inv(A).
//
            if (j+jb<=n) {
            }
            IF( J+JB.LE.N )
     $          CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                      N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                      WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
        }
    }
//
//  Apply column interchanges.
//
    for (IndexType j=n-1; j>=1; --j) {
        const IndexType jp = piv(j);
        if (jp!=j) {
            blas::swap(A(_,j), A(_,jp));
        }
    }

    work(1) = iws;
    return info;
}

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri_generic(TrMatrix<MA> &A)
{
}

//== interface for native lapack ===============================================

#ifdef CHECK_CXXLAPACK
//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri_native(GeMatrix<MA> &A, DenseVector<VP> &piv, DenseVector<VWORK> &work)
{
}

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri_native(TrMatrix<MA> &A)
{
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

//-- (ge)tri
template <typename MA, typename VP, typename VWORK>
typename GeMatrix<MA>::IndexType
tri(GeMatrix<MA> &A, DenseVector<VP> &piv, DenseVector<VWORK> &work)
{
}

//-- (tr)tri
template <typename MA>
typename GeMatrix<MA>::IndexType
tri(TrMatrix<MA> &A)
{
}

//-- forwarding ----------------------------------------------------------------

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GESV_TRI_TCC
