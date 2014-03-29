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

#ifndef FLENS_BLAS_CLOSURES_R_TCC
#define FLENS_BLAS_CLOSURES_R_TCC 1

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

//-- GeneralMatrix -------------------------------------------------------------
//
//  A += alpha*x*transpose(y)
//
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<DefaultEval<MatrixClosureRU<VX, VY> >::value
                 && IsVector<VX>::value
                 && IsVector<VY>::value
                 && IsGeneralMatrix<MA>::value,
         void>::Type
axpy(Transpose DEBUG_VAR(trans), const ALPHA &alpha,
     const MatrixClosureRU<VX, VY> &xyt, MA &A)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans==NoTrans);

    typedef typename PruneScaling<VX>::Remainder    RVX;
    typedef typename PruneScaling<VX>::ScalingType  SVX;

    const RVX &x_     = PruneScaling<VX>::getRemainder(xyt.left());
    const SVX &alpha_ = PruneScaling<VX>::getFactor(xyt.left());
    const auto &y     = xyt.right().right();

    blas::ru(alpha*alpha_, x_, y, A.impl());
}

//
//  A += alpha*x*conjTrans(y)
//
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<DefaultEval<MatrixClosureRC<VX, VY> >::value
                 && IsVector<VX>::value
                 && IsVector<VY>::value
                 && IsGeneralMatrix<MA>::value,
         void>::Type
axpy(Transpose DEBUG_VAR(trans), const ALPHA &alpha,
     const MatrixClosureRC<VX, VY> &xyh, MA &A)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans==NoTrans);

    typedef typename PruneScaling<VX>::Remainder    RVX;
    typedef typename PruneScaling<VX>::ScalingType  SVX;

    const RVX &x_     = PruneScaling<VX>::getRemainder(xyh.left());
    const SVX &alpha_ = PruneScaling<VX>::getFactor(xyh.left());
    const auto &y     = xyh.right().right().right();

    blas::rc(alpha*alpha_, x_, y, A.impl());
}

//-- SymmetricMatrix -----------------------------------------------------------
//
//  A += alpha*x*transpose(x)
//
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<DefaultEval<MatrixClosureRU<VX, VY> >::value
                 && IsVector<VX>::value
                 && IsVector<VY>::value
                 && IsSymmetricMatrix<MA>::value,
         void>::Type
axpy(Transpose DEBUG_VAR(trans), const ALPHA &alpha,
     const MatrixClosureRU<VX, VY> &xxt, MA &A)
{
    ASSERT(trans==Trans || trans==NoTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<VX>::Remainder    RVX;
    typedef typename PruneScaling<VX>::ScalingType  SVX;

    const RVX &x_     = PruneScaling<VX>::getRemainder(xxt.left());
    const SVX &alpha_ = PruneScaling<VX>::getFactor(xxt.left());

#   ifndef NDEBUG
    const auto &y     = xxt.right().right();
    ASSERT(identical(x_, y));
#   endif

    blas::r(alpha*alpha_, x_, A.impl());
}

//-- HermitianMatrix -----------------------------------------------------------
//
//  A += alpha*x*conjTrans(x)
//
template <typename ALPHA, typename VX, typename VY, typename MA>
typename RestrictTo<DefaultEval<MatrixClosureRC<VX, VY> >::value
                 && IsVector<VX>::value
                 && IsVector<VY>::value
                 && IsHermitianMatrix<MA>::value,
         void>::Type
axpy(Transpose DEBUG_VAR(trans), const ALPHA &alpha,
     const MatrixClosureRC<VX, VY> &xxh, MA &A)
{
    // Lehn: keep it simple for the moment
    ASSERT(trans==NoTrans);

    using namespace DEBUGCLOSURE;

    typedef typename PruneScaling<VX>::Remainder    RVX;
    typedef typename PruneScaling<VX>::ScalingType  SVX;

    const RVX &x_     = PruneScaling<VX>::getRemainder(xxh.left());
    const SVX &alpha_ = PruneScaling<VX>::getFactor(xxh.left());

#   ifndef NDEBUG
    const auto &y     = xxh.right().right().right();
    ASSERT(identical(x_, y));
#   endif

    blas::r(alpha*alpha_, x_, A.impl());
}


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_R_TCC
