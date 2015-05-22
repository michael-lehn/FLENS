#ifndef LU_CHECK_LU_H
#define LU_CHECK_LU_H 1

#include <flens/flens.h>
#include <flens/examples/lu/apply_perm.h>

namespace flens {

//
//  This routine is just intended to check results.  So we don't really care
//  about performance.  We just want it easy.
//
template <typename MA, typename MLU, typename VP>
double
check_LU(GeMatrix<MA> &A, const GeMatrix<MLU> &LU, const DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    const Underscore<IndexType> _;

    if (m==n) {
        typename GeMatrix<MA>::NoView L(n, n);
        L = LU.lowerUnit();
        L = L*LU.upper();
        apply_perm(p, L);
        A -= L;
    } else if (m>n) {
        typename GeMatrix<MA>::NoView LU_(m, n);
        typename GeMatrix<MA>::NoView L(m, n);
        typename GeMatrix<MA>::NoView U(n, n);
        L = LU.lowerUnit();
        U = LU(_(1,n),_(1,n)).upper();

        LU_ = L*U;
        apply_perm(p, LU_);
        A -= LU_;
    } else {
        typename GeMatrix<MA>::NoView LU_(m, n);
        typename GeMatrix<MA>::NoView L(m, m);
        typename GeMatrix<MA>::NoView U(m, n);
        L = LU(_(1,m),_(1,m)).lowerUnit();
        U = LU.upper();

        LU_ = L*U;
        apply_perm(p, LU_);
        A -= LU_;
    }


    return blas::asum(A);
}


} // namespace flens

#endif // LU_CHECK_H
