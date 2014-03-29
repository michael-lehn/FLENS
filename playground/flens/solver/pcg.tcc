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
 * Yousef Saad - Iterative methods for sparse linear systems  (2nd edition)
 * Algorithm 9.1
 *
 */

#ifndef PLAYGROUND_FLENS_SOLVER_PCG_TCC
#define PLAYGROUND_FLENS_SOLVER_PCG_TCC 1

#include <cmath>
#include <playground/flens/solver/pcg.h>

namespace flens { namespace solver {

template <typename MP, typename MA, typename VX, typename VB>
    typename RestrictTo<IsMatrix<MP>::value
                     && IsSymmetricMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
pcg(const MP &P, const MA &A, VX &&x, const VB &b,
       typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
       typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;

    Vector      Ap, r, p, z;
    ElementType alpha, beta, pNormSquare, rz, rzPrev;

    r  = b - A*x;
    z  = P*r;
    p  = z;
    rz = r*z;

    for (IndexType k=1; k<=maxIterations; k++) {

        pNormSquare = p*p;

        if (abs(pNormSquare)<=tol) {
            return 0;
        }

        Ap     = A*p;
        alpha  = rz/(p*Ap);
        x      = x + alpha*p;

        r      = r - alpha*Ap;
        z      = P*r;

        rzPrev = rz;
        rz     = r*z;
        beta   = rz/rzPrev;
        p      = beta*p + z;
    }
    return maxIterations;

}


} }// namespace solvers, flens

#endif // PLAYGROUND_FLENS_SOLVER_PCG_TCC
