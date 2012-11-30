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

#ifndef PLAYGROUND_FLENS_SOLVER_BISTABCG_TCC
#define PLAYGROUND_FLENS_SOLVER_BISTABCG_TCC 1

#include <cmath>

namespace flens { namespace solvers {

template <typename MA, typename VX, typename VB>
typename RestrictTo<IsMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename VX::IndexType>::Type
bicgstab(const MA &A, VX &x, const VB &b,
         typename ComplexTrait<typename VX::ElementType>::PrimitiveType tol,
         typename VX::IndexType maxIterations)
{
    using std::abs;

    typedef typename VX::ElementType   ElementType;
    typedef typename VX::IndexType     IndexType;
    
    VX Ap, As, r, rs, s, p;
    ElementType alpha, beta, rNormSquare, omega;
    
    const ElementType One(1);
    
    r  = b - A*x;
    rs = r;
    p  = r;
    
    rNormSquare = r*r;
    
    for (IndexType k=1; k<=maxIterations; k++) {
      
        if (abs(rNormSquare)<=tol) {
            return 0;
        }
        
        Ap    = A*p;
        alpha = blas::dotu(r, rs)/blas::dotu(Ap, rs);
        s     = r - alpha*Ap;
        As    = A*s;
        omega = blas::dotu(As, s)/blas::dotu(As, As);
        x     = x + alpha*p + omega*s;
        beta  = One/blas::dotu(r, rs);
        r     = s - omega*As;
        beta  = (beta*alpha*blas::dotu(r, rs))/omega;
        p     = beta*p - beta*omega*Ap + r;
	
        rNormSquare = blas::dotu(r, r);
    }
    return maxIterations;
}

} }// namespace solvers, flens

#endif // PLAYGROUND_FLENS_SOLVER_BISTABCG_TCC