#ifndef FLENS_EXAMPLES_CG_H
#define FLENS_EXAMPLES_CG_H 1

#include <flens/flens.h>
#include <limits>

template <typename MA, typename VX, typename VB>
int
cg(const MA &A, const VB &b, VX &&x,
   double tol = std::numeric_limits<double>::epsilon(),
   int    maxIterations = std::numeric_limits<int>::max())
{
    typedef typename VB::ElementType  ElementType;
    typedef typename VB::IndexType    IndexType;
    typedef typename VB::NoView       VectorType;

    ElementType  alpha, beta, rNormSquare, rNormSquarePrev;
    VectorType   Ap, r, p;

    r = b - A*x;
    p = r;
    rNormSquare = r*r;
    for (int k=1; k<=maxIterations; ++k) {
        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rNormSquare/(p * Ap);
        x += alpha*p;
        r -= alpha*Ap;

        rNormSquarePrev = rNormSquare;
        rNormSquare = r*r;
        beta = rNormSquare/rNormSquarePrev;
        p = beta*p + r;
    }
    return maxIterations;
}

#endif // FLENS_EXAMPLES_CG_H
