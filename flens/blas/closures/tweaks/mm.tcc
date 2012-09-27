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

#ifndef FLENS_BLAS_CLOSURES_TWEAKS_MM_TCC
#define FLENS_BLAS_CLOSURES_TWEAKS_MM_TCC 1

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
//  C = beta*D + alpha*op(A)*op(B)
//
template <typename MD, typename MA, typename MB, typename MC>
typename RestrictTo<DefaultEval<MatrixClosureMM<MD, MA, MB> >::value
                 && IsMatrix<MD>::value
                 && IsMatrix<MA>::value
                 && IsMatrix<MB>::value
                 && IsGeneralMatrix<MC>::value,
         void>::Type
copy(Transpose trans, const MatrixClosureMM<MD, MA, MB> &DpAB, MC &C)
{
    FLENS_BLASLOG_BEGIN_COPY(DpAB, C);

    using namespace DEBUGCLOSURE;
//
//  check if C form rhs and lhs are identical
//
    typedef typename PruneScaling<MD>::Remainder    RMD;

    const RMD &D = PruneScaling<MD>::getRemainder(DpAB.left());

    typedef typename MC::Impl::ElementType  TC;
    const TC  One(1);

    if (trans!=NoTrans || !identical(D, C.impl())) {
        copySum(trans, DpAB.left(), One, DpAB.right(), C.impl());
        FLENS_BLASLOG_END;
        return;
    }
//
//  get factor  beta
//
    typedef typename PruneScaling<MD>::ScalingType  SMD;
    const SMD &beta = PruneScaling<MD>::getFactor(DpAB.left());
//
//  Rest gets done by the mm switch
//
    const auto &A = DpAB.right().left();
    const auto &B = DpAB.right().right();

    if (trans==NoTrans || trans==Conj) {
        mmSwitch(trans, trans, One, A, B, beta, C.impl());
    } else {
        mmSwitch(trans, trans, One, B, A, beta, C.impl());
    }

    FLENS_BLASLOG_END;
}


} } // namespace blas, flens

#endif // FLENS_BLAS_CLOSURES_TWEAKS_MM_TCC

