/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_BLAS_LEVEL3_SM_TCC
#define FLENS_BLAS_LEVEL3_SM_TCC

#include <flens/blas/closures/closures.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

namespace flens { namespace blas {

//-- trsm
template <typename ALPHA, typename MA, typename MB>
typename RestrictTo<IsTrMatrix<MA>::value
                 && IsGeMatrix<MB>::value,
             void>::Type
sm(Side             side,
   Transpose        transA,
   const ALPHA      &alpha,
   const MA         &A,
   MB               &&B)
{
    typedef typename RemoveRef<MA>::Type MatrixA;
    typedef typename RemoveRef<MB>::Type MatrixB;

#   ifndef NDEBUG
    ASSERT(B.order()==A.order());
    if (side==Left) {
        assert(A.dim()==B.numRows());
    } else {
        assert(B.numCols()==A.dim());
    }
#   endif

#   ifdef HAVE_CXXBLAS_TRSM
    cxxblas::trsm(A.order(), side, A.upLo(),
                  transA, A.diag(),
                  B.numRows(), B.numCols(),
                  alpha,
                  A.data(), A.leadingDimension(),
                  B.data(), B.leadingDimension());
#   else
    ASSERT(0);
#   endif
}

} } // namespace blas, flens

#endif // FLENS_BLAS_LEVEL3_SM_TCC
