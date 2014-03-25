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
 * Wenwu Chen and Bill Poirer 
 * - Parallel implementation of efficient preconditioned linear solver for 
 *   grid-based applications in chemical physics. II: QMR linear solver
 * Journal of Computational Physics 219 (2006) 198–209
 *
 */

#ifndef PLAYGROUND_FLENS_SOLVER_TFQMR_TCC
#define PLAYGROUND_FLENS_SOLVER_TFQMR_TCC 1

#include <cmath>

namespace flens { namespace solver {

template <typename MA, typename VX, typename VB>
    typename RestrictTo<IsGeneralMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
tfqmr(const MA &A, VX &&x, const VB &b,
      typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
      typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;
   

    Vector s, d, v, p, Atp, q, r;
    
    ElementType beta, gamma, delta;
    ElementType epsilon, eta, theta, rho;
    ElementType gammaold, thetaold, rhoold;

    const ElementType One(1);

    gamma = One;
    theta = One;

    r = b - A*x;
    v = transpose(A)*r;
    r = v;
    
    rho = blas::nrm2(v);
    
    for(IndexType i = 0; i < maxIterations; ++i)
    {

        if( abs(r*r)<tol ) {
            return i;
            
        }
        
        ASSERT( abs(rho)!=0 );

        v = (One/rho)*v;

        if ( i==0 ) {

            delta = (One/rho)*(v*b);
        } else {
        
            delta = v*v;
            
        }
        ASSERT( abs(delta)!=0 );

        if(i == 0) {
            q = (One/rho)*b;
        } else {
            q = (-rho*delta/epsilon)*q + v;
        }

        p       = A*q;
        Atp     = transpose(A)*p;
        epsilon = q*Atp;

        ASSERT( abs(epsilon)!=0 );

        beta = epsilon/delta;
        v    = -beta*v + Atp;

        rhoold = rho;
        rho    = blas::nrm2(v);

        thetaold = theta;
        theta    = rho/(gamma*abs(beta));

        gammaold = gamma;
        gamma    = One/(sqrt(One+theta*theta));

        ASSERT( abs(gamma)!=0 );

        if(i == 0)
        {
            eta = rhoold*gamma*gamma/beta;

            d = eta*q;
            s = eta*Atp;
        }
        else
        {

            ElementType kappa = pow(gamma/gammaold,2);
            eta = -eta*rhoold*kappa/beta;
  
            ElementType lambda = pow(thetaold*gamma,2);          
            d = lambda*d + eta*q;
            s = lambda*s + eta*Atp;
        }

        x +=  d;
        r -=  s;
    }
    return maxIterations;
}


template <typename MA, typename VX, typename VB>
    typename RestrictTo<IsSymmetricMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
tfqmr(const MA &A, VX &&x, const VB &b,
      typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
      typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;
   

    Vector s, d, v, p, q, r;
    
    ElementType beta, gamma, delta;
    ElementType epsilon, eta, theta, rho;
    ElementType gammaold, thetaold, rhoold;

    const ElementType One(1);

    gamma = One;
    theta = One;

    r = b - A*x;
    v = r;

    rho = blas::nrm2(v);
    
    for(IndexType i = 0; i < maxIterations; ++i)
    {

        if( abs(r*r)<tol ) {
            return i;
            
        }
                 ASSERT( abs(rho)!=0 );

        v = (One/rho)*v;

        if ( i==0 ) {

            delta = (One/rho)*(v*b);
        } else {
        
            delta = v*v;
            
        }
        ASSERT( abs(delta)!=0 );

        if(i == 0) {
            q = (One/rho)*b;
        } else {
            q = (-rho*delta/epsilon)*q + v;
        }

        p       = A*q;
        epsilon = q*p;

        ASSERT( abs(epsilon)!=0 );

        beta = epsilon/delta;
        v    = -beta*v + p;

        rhoold = rho;
        rho    = blas::nrm2(v);

        thetaold = theta;
        theta    = rho/(gamma*abs(beta));

        gammaold = gamma;
        gamma    = One/(sqrt(One+theta*theta));

        ASSERT( abs(gamma)!=0 );

        if(i == 0)
        {
            eta = rhoold*gamma*gamma/beta;

            d = eta*q;
            s = eta*p;
        }
        else
        {

            ElementType kappa = pow(gamma/gammaold,2);
            eta = -eta*rhoold*kappa/beta;
  
            ElementType lambda = pow(thetaold*gamma,2);          
            d = lambda*d + eta*q;
            s = lambda*s + eta*p;
        }

        x +=  d;
        r -=  s;

    }
    return maxIterations;
}


template <typename MP, typename MA, typename VX, typename VB>
    typename RestrictTo<IsGeneralMatrix<MP>::value
                     && IsSymmetricMatrix<MA>::value
                     && IsDenseVector<VX>::value
                     && IsDenseVector<VB>::value,
             typename RemoveRef<VX>::Type::IndexType>::Type
ptfqmr(const MP &P, const MA &A, VX &&x, const VB &b,
      typename ComplexTrait<typename RemoveRef<VX>::Type::ElementType>::PrimitiveType tol,
      typename RemoveRef<VX>::Type::IndexType maxIterations)
{
    using std::abs;

    typedef typename RemoveRef<VX>::Type   VectorX;
    typedef typename VectorX::NoView       Vector;
    typedef typename VectorX::IndexType    IndexType;
    typedef typename VectorX::ElementType  ElementType;
   

    Vector s, d, v, z, p, q, r;
    
    ElementType alpha, beta, gamma, delta;
    ElementType epsilon, eta, theta, xi, rho;
    ElementType gammaold, thetaold, rhoold;

    const ElementType One(1);

    alpha = One;
    gamma = One;
    theta = One;

    z = P*b;
    r = b - A*x;
    v = r;

    rho = blas::nrm2(v);
    xi  = blas::nrm2(v);

    
    for(IndexType i = 0; i < maxIterations; ++i)
    {

        if( abs(r*r)<tol ) {
            return i;
            
        }
        
        ASSERT( abs(rho)!=0 );
        ASSERT( abs(xi)!=0 );
        
        alpha = alpha*(xi/rho);

        v = (One/rho)*v;
        z = (One/xi)*z;

        delta = v*z;

        ASSERT( abs(delta)!=0 );

        if(i == 0) {
            q = z;
        } else {
            q = (-rho*delta/epsilon)*q + z;
        }

        p       = alpha*A*q;
        epsilon = q*p;

        ASSERT( abs(epsilon)!=0 );

        beta = epsilon/delta;
        v    = -beta*v + p;

        rhoold = rho;
        rho    = blas::nrm2(v);

        z  = (One/alpha)*P*v;
        
        xi = blas::nrm2(z);

        thetaold = theta;
        theta    = rho/(gamma*abs(beta));

        gammaold = gamma;
        gamma    = One/(sqrt(One+theta*theta));

        ASSERT( abs(gamma)!=0 );

        if(i == 0)
        {
            eta = rhoold*gamma*gamma/beta;

            d = eta*alpha*q;
            s = eta*p;
        }
        else
        {
            ElementType lambda = pow(thetaold*gamma, 2);
            ElementType kappa  = pow(gamma/gammaold, 2);
            eta = -eta*rhoold*kappa/beta;
            
            d = lambda*d + eta*alpha*q;
            s = lambda*s + eta*p;
        }

        x +=  d;
        r -=  s;
    }
    return maxIterations;
}

} }// namespace solver, flens

#endif // PLAYGROUND_FLENS_SOLVER_TFQMR_TCC