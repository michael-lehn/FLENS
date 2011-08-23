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
       SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
 *
 *  -- LAPACK auxiliary routine (version 3.3.1)                        --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2009                                                      --
 */

#ifndef FLENS_LAPACK_EIG_LAHR2_TCC
#define FLENS_LAPACK_EIG_LAHR2_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//-- forwarding ----------------------------------------------------------------
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
void
lahr2(IndexType k, IndexType nb, MA &&A, VTAU &&tau, MTR &&Tr, MY &&Y)
{
    lahr2(k, nb, A, tau, Tr, Y);
}

//-- lahr2 ---------------------------------------------------------------------
template <typename IndexType, typename MA, typename VTAU,
          typename MTR, typename MY>
void
lahr2(IndexType k, IndexType nb, GeMatrix<MA> &A, DenseVector<VTAU> &tau,
      GeMatrix<MTR> &_Tr, GeMatrix<MY> &Y)
{
    // TODO: assertions for ranges of k and nb
    ASSERT(k<A.numRows());
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numCols()==A.numRows()-k+1);
    ASSERT(tau.length()==nb);
    ASSERT(_Tr.numRows()==_Tr.numCols());
    ASSERT(_Tr.numRows()==nb);
    ASSERT(Y.numRows()==A.numRows());
    ASSERT(Y.numCols()==nb);

    typedef typename GeMatrix<MA>::ElementType T;

    const T         One(1);
    const IndexType n = A.numRows();

}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_LAHR2_TCC
