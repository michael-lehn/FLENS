/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_BLAS_CLOSURES_COPY_TCC
#define FLENS_BLAS_CLOSURES_COPY_TCC 1

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

//------------------------------------------------------------------------------
//
//  y = beta*z + alpha*op(A)*x
//
template <typename VZ, typename MA, typename VX, typename VY>
typename RestrictTo<DefaultEval<VectorClosureMV<VZ, MA, VX> >::value
                 && IsVector<VZ>::value
                 && IsMatrix<MA>::value
                 && IsVector<VX>::value,
         void>::Type
copy(const VectorClosureMV<VZ, MA, VX> &zpAx, Vector<VY> &y)
{
    FLENS_BLASLOG_BEGIN_COPY(zpAx, y);

    using namespace DEBUGCLOSURE;
//
//  check if z form rhs and y from lhs are identical
//
    typedef typename PruneScaling<VZ>::Remainder    RVZ;

    const RVZ &z = PruneScaling<VZ>::getRemainder(zpAx.left());

    if (!identical(z, y.impl())) {
        typedef typename VY::Impl::ElementType  TY;
        const TY  One(1);
//
//      z and y are not identical.  Compute:
//        (1) y  = beta*z
//        (2) y += alpha*op(A)*x
//
        copySum(zpAx.left(), One, zpAx.right(), y.impl());
        FLENS_BLASLOG_END;
        return;
    }
//
//  z and y are identical.  Get factor beta ...
//
    typedef typename PruneScaling<VZ>::ScalingType  SVZ;
    const SVZ &beta = PruneScaling<VZ>::getFactor(zpAx.left());
//
//  ... and let the mv switch do the rest.
//
    typedef typename MA::Impl::ElementType  TA;
    const auto &A = zpAx.right().left();
    const auto &x = zpAx.right().right();

    mvSwitch(NoTrans, TA(1), A, x, beta, y.impl());

    FLENS_BLASLOG_END;
}
} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_COPY_TCC
