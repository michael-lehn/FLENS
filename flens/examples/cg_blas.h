#ifndef FLENS_EXAMPLES_CG_BLAS_H
#define FLENS_EXAMPLES_CG_BLAS_H 1

#include <flens/flens.h>
#include <cxxstd/limits.h>

template <typename MA, typename VX, typename VB>
int
cg(const MA &A, const VB &b, VX &&x,
   double tol = std::numeric_limits<double>::epsilon(),
   int    maxIterations = std::numeric_limits<int>::max())
{
    using namespace flens;

    typedef typename VB::ElementType  ElementType;
    typedef typename VB::IndexType    IndexType;
    typedef typename VB::NoView       VectorType;

    ElementType  alpha, beta, rNormSquare, rNormSquarePrev;
    VectorType   Ap, r, p;

    const ElementType  Zero(0), One(1);

///
/// `r = b - A*x;`
///
    blas::copy(b, r);
    blas::mv(NoTrans, -One, A, x, One, r);

///
/// `p = r;`
///
    blas::copy(r, p);

///
/// `rNormSquare = r*r;`
///
    rNormSquare = blas::dot(r, r);

    for (int k=1; k<=maxIterations; ++k) {
        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }

///
///     `Ap = A*p;`
///
        blas::mv(NoTrans, One, A, p, Zero, Ap);

///
///     `alpha = rNormSquare/(p * Ap);`
///
        alpha = rNormSquare/blas::dot(p, Ap);

///
///     `x += alpha*p;`
///
        blas::axpy(alpha, p, x);

///
///     `r -= alpha*Ap;`
///
        blas::axpy(-alpha, Ap, r);

        rNormSquarePrev = rNormSquare;

///
///     `rNormSquare = r*r;`
///
        rNormSquare = blas::dot(r, r);

        beta = rNormSquare/rNormSquarePrev;

///
///     `p = beta*p + r;`
///
        blas::scal(beta, p);
        blas::axpy(One, r, p);
    }
    return maxIterations;
}

#endif // FLENS_EXAMPLES_CG_BLAS_H
