#ifndef LU_LU_TILED_H
#define LU_LU_TILED_H 1

#ifndef LU_BS
#define LU_BS 96
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_blk_with_operators.h>
#include <flens/examples/lu/tile.h>

namespace flens {

template <typename MA, typename VP>
typename DenseVector<VP>::IndexType
lu_tiled(TiledCopy<MA> &A, DenseVector<VP> &p)
{
    typedef typename TiledCopy<MA>::ElementType  T;
    typedef typename TiledCopy<MA>::IndexType    IndexType;

    const IndexType  m  = A.numTileRows();
    const IndexType  n  = A.numTileCols();
    const IndexType  bs = A.blockSize();
    const IndexType  mn = std::min(m, n);

    const T One(1);

    const Underscore<IndexType>  _;

    IndexType info = 0;

    for (IndexType i=1; i<=mn; ++i) {
        auto A_ii = A(i,i);
        auto p_   = p(_(bs*(i-1)+1,bs*(i-1)+A_ii.numRows()));

        info = lu_blk(A_ii, p_);

        if (info) {
            return info+(i-1)*bs;
        }

        for (IndexType j=1; j<=n; ++j) {
            if (i!=j) {
                auto A_ij = A(i,j);
                apply_perm_inv(p_, A_ij);
            }
        }
        p_ += bs*(i-1);

        for (IndexType k=i+1; k<=m; ++k) {
            blas::sm(Right, NoTrans, One, A_ii.upper(), A(k,i));
        }
        for (IndexType l=i+1; l<=n; ++l) {
            blas::sm(Left, NoTrans, One, A_ii.lowerUnit(), A(i,l));
        }
        for (IndexType k=i+1; k<=m; ++k) {
            for (IndexType l=i+1; l<=n; ++l) {
                A(k,l) -= A(k,i)*A(i,l);
            }
        }
    }
    return 0;
}

} // namespace flens

#endif // LU_LU_TILED_H
