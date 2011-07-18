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

#ifndef EXTERNAL_ITERATIVESOLVERS_PCG_TCC
#define EXTERNAL_ITERATIVESOLVERS_PCG_TCC 1

#include <external/iterativesolvers/cg.h>

namespace flens {

template <typename Prec, typename MA, typename VA, typename VB, typename TOL,
          typename IndexType>
IndexType
pcg(const Prec &B, const MA &A, VA &x, const VB &b, TOL tol,
    IndexType maxIterations)
{
    typename cg_<VB>::AuxVector  r, q, p, Ap;
    typename cg_<VB>::T          pNormSquare, alpha, beta, rq, rqPrev;

    r = b - A*x;
    q = B*r;
    p = q;
    rq = r*q;
    for (IndexType k=1; k<=maxIterations; k++) {
        pNormSquare = p*p;
        if (pNormSquare<=tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rq/(p*Ap);
        x += alpha*p;

        r -= alpha*Ap;
        q = B*r;

        rqPrev = rq;
        rq = r*q;
        beta = rq/rqPrev;
            p = beta*p + q;
    }
    return maxIterations;
}

} // namespace flens

#endif // EXTERNAL_ITERATIVESOLVERS_PCG_TCC
