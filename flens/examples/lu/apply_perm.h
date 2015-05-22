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
    const IndexType n = A.numCols();

    const IndexType bs = 32;
    const IndexType nbs = (n/bs)*bs;

    // reverse row interchanges
    if (nbs!=0) {
        for (IndexType j=1; j<=nbs; j+=bs) {
            for (IndexType i=m; i>=1; --i) {
                if (i!=p(i)) {
                    blas::swap(A(i,_(j,j+bs-1)), A(p(i),_(j,j+bs-1)));
                }
            }
        }
    }
    if (nbs!=n) {
        for (IndexType j=1; j<=nbs; j+=bs) {
            for (IndexType i=m; i>=1; --i) {
                if (i!=p(i)) {
                    blas::swap(A(i,_(nbs+1,n)), A(p(i),_(nbs+1,n)));
                }
            }
        }
    }
    /*
    for (IndexType i=m; i>=1; --i) {
        if (i!=p(i)) {
            blas::swap(A(i,_), A(p(i),_));
        }
    }
    */
}

} // namespace flens

#endif // LU_APPLY_PERM_H
