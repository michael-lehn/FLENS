#ifndef LU_LU_UNBLK_H
#define LU_LU_UNBLK_H 1

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/limits.h>
#include <flens/flens.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_unblk(GeMatrix<MA> &A, DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = std::min(m, n);

    const T Zero(0);
    const T One(1);

    const Underscore<IndexType>  _;

    // Compute threshold for stable scaling
    const T eps   = std::numeric_limits<T>::epsilon() * T(0.5);
    T safeMin     = std::numeric_limits<T>::min();
    const T small = One / std::numeric_limits<T>::max();
    if (small>=safeMin) {
        safeMin = small*(One+eps);
    }

    for (IndexType i=1; i<=mn; ++i) {
///
///     The partitioning is realized by referencing matrix blocks.  Referencing
///     means we either have matrix/vector views (i.e.`A_1`, `a12`, `a21`,
///     `A22`) or a C++ reference for scalars (i.e. `a11`).  So changing any
///     of these changes the original matrix.
///
        // Partition sub-matrix A(i:m,i:n)
        const auto  A_1 = A(_( i, m),        i);  // first column of sub_matrix
        const auto &a11 = A(       i,        i);
        const auto  a12 = A(       i, _(i+1,n));
        auto        a21 = A(_(i+1,m),        i);
        auto        A22 = A(_(i+1,m), _(i+1,n));

        // Find pivot index
        IndexType l = blas::iamax(A_1) +i- 1;
        p(i) = l;

        // Interchange rows of A
        if (p(i)!=i) {
            blas::swap(A(i,_), A(l,_));
        }

        // Only abort if diagonal element is *exactly* zero
        if (a11==Zero) {
            return i;
        }

        // Choose stable scaling
        if (std::abs(a11)<safeMin) {
            a21 /= a11;
        } else {
            a21 *= One/a11;
        }

        // Rank 1 update
        A22 -= a21*transpose(a12);
    }
    return 0;
}

} // namespace flens

#endif // LU_LU_UNBLK_H
