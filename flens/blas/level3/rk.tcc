/*
 *   Copyright (c) 2010,2015 Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_RK_TCC
#define FLENS_BLAS_LEVEL3_RK_TCC

#include <flens/blas/closures/closures.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>
#include <ulmblas/cxxblas.h>

namespace flens { namespace blas {

//-- herk
template <typename ALPHA, typename MA, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsHeMatrix<MC>::value,
         void>::Type
rk(Transpose        transposeA,
   const ALPHA      &alpha,
   const MA         &A,
   const BETA       &beta,
   MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType IndexType;

    const bool transA = (transposeA==Trans || transposeA==ConjTrans);
    const bool lowerC = (C.upLo()==Lower);

    const IndexType n = (!transA) ? A.numRows() : A.numCols();
    const IndexType k = (!transA) ? A.numCols() : A.numRows();

    ASSERT(transposeA!=Conj && transposeA!=Trans);
    ASSERT((beta==BETA(0)) || (C.dim()==n));
    if (C.dim()!=n) {
        C.resize(n, n);
    }

    ASSERT(C.dim()==((!transA) ? A.numRows() : A.numCols()));

    cxxblas::herk(n, k, alpha,
                  transA, A.data(), A.strideRow(), A.strideCol(),
                  beta,
                  lowerC, C.data(), C.strideRow(), C.strideCol());
}

//-- syrk
template <typename ALPHA, typename MA, typename BETA, typename MC>
typename RestrictTo<IsGeMatrix<MA>::value
                 && IsSyMatrix<MC>::value,
         void>::Type
rk(Transpose        transposeA,
   const ALPHA      &alpha,
   const MA         &A,
   const BETA       &beta,
   MC               &&C)
{
    typedef typename RemoveRef<MC>::Type MatrixC;
    typedef typename MatrixC::IndexType IndexType;

    const bool transA = (transposeA==Trans || transposeA==ConjTrans);
    const bool lowerC = (C.upLo()==Lower);

    const IndexType n = (!transA) ? A.numRows() : A.numCols();
    const IndexType k = (!transA) ? A.numCols() : A.numRows();

    //ASSERT(transposeA!=Conj && transposeA!=ConjTrans);
    ASSERT((beta==BETA(0)) || (C.dim()==n));
    if (C.dim()!=n) {
        C.resize(n, n);
    }

    ASSERT(C.dim()==((!transA) ? A.numRows() : A.numCols()));

    cxxblas::syrk(n, k, alpha,
                  transA, A.data(), A.strideRow(), A.strideCol(),
                  beta,
                  lowerC, C.data(), C.strideRow(), C.strideCol());
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_RK_TCC
