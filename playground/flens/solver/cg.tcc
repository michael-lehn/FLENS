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
 * Algorithm 6.18
 *
 */

#ifndef PLAYGROUND_FLENS_SOLVER_CG_TCC
#define PLAYGROUND_FLENS_SOLVER_CG_TCC 1

#include <cmath>
#include <playground/flens/solver/cg.h>

namespace flens { namespace solver {

template <typename MA, typename VX, typename VB>
    typename RestrictTo<IsSymmetricMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
cg(const MA &A, VX &&x, const VB &b,
    typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
    typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;
    
    Vector Ap, r, p;
    ElementType alpha, beta, rNormSquare, rNormSquarePrev;
    
    r = b - A*x;
    p = r;
    rNormSquare = r*r;
    for (IndexType k=1; k<=maxIterations; k++) {

        if (abs(rNormSquare)<=tol) {
            return 0;
        }
        
        Ap    = A*p;
        alpha = rNormSquare/(p*Ap);
        x     = x + alpha*p;
        r     = r - alpha*Ap;

        rNormSquarePrev = rNormSquare;
        rNormSquare     = r*r;
        
        beta = rNormSquare/rNormSquarePrev;
        p    = beta*p + r;
    }
    return maxIterations;
}

template <typename MA, typename VX, typename VB>
    typename RestrictTo<IsGeneralMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
cg(const MA &A, VX &&x, const VB &b,
    typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
    typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;
   
    Vector      Ap, AtAp, r, p;
    ElementType alpha, beta, rNormSquare, rNormSquarePrev;
    
    p = b - A*x;
    r = transpose(A)*p;
    p = r;
    
    rNormSquare = r*r;
    for (IndexType k=1; k<=maxIterations; k++) {
      
        if ( abs(rNormSquare)<=tol ) {
            return 0;
        }

        Ap    = A*p;
        AtAp  = transpose(A)*Ap;
        alpha = rNormSquare/(p*AtAp);
        x     = x + alpha*p;
        r     = r - alpha*AtAp;

        rNormSquarePrev = rNormSquare;
        rNormSquare     = r*r;

        beta  = rNormSquare/rNormSquarePrev;
        p     = beta*p + r;
    }
    return maxIterations;
}

} }// namespace solver, flens

#endif // PLAYGROUND_FLENS_SOLVER_CG_TCC
