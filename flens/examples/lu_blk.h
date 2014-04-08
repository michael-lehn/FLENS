#include <algorithm>
#include <flens/flens.h>
#include <flens/examples/lu_unblk.h>

#ifndef BS
#define BS 32
#endif


namespace flens {

template <typename MA>
typename GeMatrix<MA>::IndexType
lu_blk(GeMatrix<MA> &A)
{
    using flens::min;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = min(m, n);

    const ElementType One(1);

    const Underscore<IndexType>  _;

    IndexType info;

    for (IndexType j=1; j<=mn; j+=BS) {
        const IndexType bs = min(BS, m-j, n-j);

        auto A11 = A(_(j,j+bs-1), _(j,j+bs-1));
        auto A12 = A(_(j,j+bs-1), _(j+bs,n));

        auto A21 = A(_(j+bs,m), _(j,j+bs-1));
        auto A22 = A(_(j+bs,m), _(j+bs,n));

        info = lu_unblk(A11);

        if (info) {
            return info+j-1;
        }

        // Compute A21*inv(A11.upper())
        blas::sm(Right, NoTrans, One, A11.upper(), A21);

        // Compute inv(A11.lowerUnit())*A12
        blas::sm(Left, NoTrans, One, A11.lowerUnit(), A12);

        A22 -= A21*A12;
    }
    return 0;
}

} // namespace flens
