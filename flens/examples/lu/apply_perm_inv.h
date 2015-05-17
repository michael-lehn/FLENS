#ifndef LU_APPLY_PERM_INV_H
#define LU_APPLY_PERM_INV_H 1

#include <flens/flens.cxx>

namespace flens {

// Compute A = P^{-1}*A = P^T*A
template <typename VP, typename MA>
void
apply_perm_inv(const DenseVector<VP> &p, GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
    const Underscore<IndexType> _;

    const IndexType m = p.length();

    // apply row interchanges
    for (IndexType i=1; i<=m; ++i) {
        if (i!=p(i)) {
            blas::swap(A(i,_), A(p(i),_));
        }
    }
}

} // namespace flens

#endif // LU_APPLY_PERM_INV_H
