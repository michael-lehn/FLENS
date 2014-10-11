#include <cxxstd/algorithm.h>
#include <flens/flens.h>

namespace flens {

template <typename MA>
typename GeMatrix<MA>::IndexType
lu_unblk(GeMatrix<MA> &A)
{
    using std::min;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = min(m, n);

    const ElementType Zero(0);

    const Underscore<IndexType>  _;

    for (IndexType j=1; j<=mn; ++j) {
        const auto a11 = A(       j,        j);
        const auto a12 = A(       j, _(j+1,n));
        auto       a21 = A(_(j+1,m),        j);
        auto       A22 = A(_(j+1,m), _(j+1,n));

        if (a11==Zero) {
            return j;
        }

        a21 /= a11;
        A22 -= a21*transpose(a12);
    }
    return 0;
}

} // namespace flens
