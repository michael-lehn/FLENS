#ifndef LU_LU_BLK_H
#define LU_LU_BLK_H 1

#ifndef LU_BS
#define LU_BS 96
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_unblk_with_operators.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_blk(GeMatrix<MA> &A, DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = std::min(m, n);

    const T One(1);

    const Underscore<IndexType>  _;

    IndexType info = 0;

    for (IndexType i=1; i<=mn; i+=LU_BS) {
        const IndexType bs = std::min(std::min(LU_BS, m-i+1), n-i+1);

        // Partitions of A
        auto A10 = A(_(   i, i+bs-1), _(   1,    i-1));
        auto A11 = A(_(   i, i+bs-1), _(   i, i+bs-1));
        auto A12 = A(_(   i, i+bs-1), _(i+bs,      n));
        auto A21 = A(_(i+bs,      m), _(   i, i+bs-1));
        auto A22 = A(_(i+bs,      m), _(i+bs,      n));

        // Part of the pivot vector for rows of A11
        auto p_  = p(_(i,i+bs-1));

        // Compute LU factorization of A11
        info = lu_unblk(A11, p_);

        if (info) {
            // All values in column info of A11 are *exactly* zero.
            return info+i- 1;
        }

        // Apply permutation to A10 and A12
        apply_perm_inv(p_, A10);
        apply_perm_inv(p_, A12);

        // Update p
        p_ += i-1;

        // Use triangular solver for A21 = A11.upper()*A21
        blas::sm(Right, NoTrans, One, A11.upper(), A21);

        // Use triangular solver for A12 = A11.lowerUnit()*A12
        blas::sm(Left, NoTrans, One, A11.lowerUnit(), A12);

        // Update A22 with matrix-product A22 = A22 - A21*A12
        A22 -= A21*A12;
    }
    return 0;
}

} // namespace flens

#endif // LU_LU_BLK_H
