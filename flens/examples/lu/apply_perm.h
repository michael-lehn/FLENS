#ifndef LU_APPLY_PERM_H
#define LU_APPLY_PERM_H 1

#include <flens/flens.cxx>

namespace flens {

// Compute A = P*A
template <typename VP, typename MA>
void
apply_perm(const DenseVector<VP> &p, GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::IndexType    IndexType;
    const Underscore<IndexType> _;

    const IndexType m = p.length();

    // reverse row interchanges
    for (IndexType i=m; i>=1; --i) {
        if (i!=p(i)) {
            blas::swap(A(i,_), A(p(i),_));
        }
    }
}

} // namespace flens

#endif // LU_APPLY_PERM_H
