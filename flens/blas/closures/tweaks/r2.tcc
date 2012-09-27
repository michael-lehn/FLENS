/*
 *   Copyright (c) 2011, Michael Lehn
 *
 *   All rights reserved.
 
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

#ifndef FLENS_BLAS_CLOSURES_R2_TCC
#define FLENS_BLAS_CLOSURES_R2_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/closures/closures.h>
#include <flens/blas/level1/level1.h>
#include <flens/blas/level2/level2.h>
#include <flens/blas/level3/level3.h>
#include <flens/typedefs.h>

#ifdef FLENS_DEBUG_CLOSURES
#   include <flens/blas/blaslogon.h>
#else
#   include <flens/blas/blaslogoff.h>
#endif

namespace flens { namespace blas {

//-- SymmetricMatrix -----------------------------------------------------------
//
//  A += a1*x1*transpose(y1) + a2*y2*transpose(x2)
//
template <typename ALPHA,
          typename VX1, typename VY1,
          typename VY2, typename VX2,
          typename MA>
typename RestrictTo<DefaultEval<MatrixClosureR2U<VX1, VY1, VY2, VX2> >::value
                 && IsVector<VX1>::value
                 && IsVector<VY1>::value
                 && IsVector<VY2>::value
                 && IsVector<VX2>::value
                 && IsSymmetricMatrix<MA>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosureR2U<VX1, VY1, VY2, VX2> &xyt_yxt, MA &A)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans==NoTrans);

    using namespace DEBUGCLOSURE;

    const auto &a1 = PruneScaling<VX1>::getFactor(xyt_yxt.left().left());
    const auto &x1 = PruneScaling<VX1>::getRemainder(xyt_yxt.left().left());
    const auto &y1 = xyt_yxt.left().right().right();

    const auto &a2 = PruneScaling<VY2>::getFactor(xyt_yxt.right().left());
    const auto &y2 = PruneScaling<VY2>::getRemainder(xyt_yxt.right().left());
    const auto &x2 = xyt_yxt.right().right().right();

    if (identical(x1,x2) && identical(y1,y2) && (a1==a2)) {
        blas::r2(alpha*a1, x1, y1, A.impl());
        return;
    } else {
        blas::axpy(trans, alpha, xyt_yxt.left(), A.impl());
        blas::axpy(trans, alpha, xyt_yxt.right(), A.impl());
    }
}

//-- HermitianMatrix -----------------------------------------------------------
//
//  A += a1*x1*conjTrans(y1) + conjugate(a2)*y2*conjTrans(x2)
//
template <typename ALPHA,
          typename VX1, typename VY1,
          typename VY2, typename VX2,
          typename MA>
typename RestrictTo<DefaultEval<MatrixClosureR2C<VX1, VY1, VY2, VX2> >::value
                 && IsVector<VX1>::value
                 && IsVector<VY1>::value
                 && IsVector<VY2>::value
                 && IsVector<VX2>::value
                 && IsHermitianMatrix<MA>::value,
         void>::Type
axpy(Transpose trans, const ALPHA &alpha,
     const MatrixClosureR2C<VX1, VY1, VY2, VX2> &xyh_yxh, MA &A)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans==NoTrans);

    using namespace DEBUGCLOSURE;

    const auto &a1 = PruneScaling<VX1>::getFactor(xyh_yxh.left().left());
    const auto &x1 = PruneScaling<VX1>::getRemainder(xyh_yxh.left().left());
    const auto &y1 = xyh_yxh.left().right().right().right();

    const auto &a2 = PruneScaling<VY2>::getFactor(xyh_yxh.right().left());
    const auto &y2 = PruneScaling<VY2>::getRemainder(xyh_yxh.right().left());
    const auto &x2 = xyh_yxh.right().right().right().right();

    if (identical(x1,x2) && identical(y1,y2) && (a1==cxxblas::conjugate(a2))) {
        blas::r2(alpha*a1, x1, y1, A.impl());
        return;
    } else {
        blas::axpy(trans, alpha, xyh_yxh.left(), A.impl());
        blas::axpy(trans, alpha, xyh_yxh.right(), A.impl());
    }
}

} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_R2_TCC
