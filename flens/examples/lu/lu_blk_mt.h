#ifndef LU_LU_BLK_MT_H
#define LU_LU_BLK_MT_H 1

#ifndef LU_BS
#define LU_BS 96
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>
#include <thread>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_unblk_with_operators.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_blk(int BlockSize, GeMatrix<MA> &A, DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = std::min(m, n);

    const T One(1);

    const Underscore<IndexType>  _;

    IndexType info = 0;

    if ((BlockSize<=1) || (BlockSize>=mn)) {
//
//      Use unblocked code.
//
        info = lu_unblk(A, p);
    } else {
        for (IndexType i=1; i<=mn; i+=BlockSize) {
            const IndexType bs = std::min(std::min(BlockSize, m-i+1), n-i+1);

            // Partitions of A
            auto A10 = A(_(   i, i+bs-1), _(   1,    i-1));
            auto A_1 = A(_(   i,      m), _(   i, i+bs-1));
            auto A11 = A(_(   i, i+bs-1), _(   i, i+bs-1));
            auto A12 = A(_(   i, i+bs-1), _(i+bs,      n));
            auto A22 = A(_(i+bs,      m), _(i+bs,      n));

            // Part of the pivot vector for rows of A11
            auto p_  = p(_(i,m));

            // Compute LU factorization of A_1
            info = lu_unblk(A_1, p_);

            if (info) {
                // All values in column info of A11 are *exactly* zero.
                return info+i-1;
            }

            // Apply permutation to A10 and A12
            apply_perm_inv(p_, A10);
            apply_perm_inv(p_, A12);

            // Update p
            p_ += i-1;

            // Use triangular solver for A12 = A11.lowerUnit()*A12
            blas::sm(Left, NoTrans, One, A11.lowerUnit(), A12);

            // Update A22 with matrix-product A22 = A22 - A21*A12
            blas::mm(NoTrans, NoTrans, -One, A21, A12, One, A22);
        }
    }
    return 0;
}


template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_blk_mt(GeMatrix<MA> &A, DenseVector<VP> &p)
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
        info = lu_blk(32, A11, p_);

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
        //                   and for A12 = A11.lowerUnit()*A12
        auto trsm0 = std::thread( [&]{ blas::sm(Right, NoTrans, One, A11.upper(), A21);    });
        auto trsm1 = std::thread( [&]{ blas::sm(Left, NoTrans, One, A11.lowerUnit(), A12); });
        trsm0.join();
        trsm1.join();

        // Update A22 with matrix-product A22 = A22 - A21*A12
        const IndexType N = A22.numCols();

        auto A22_left  = A22(_,_(    1, N/2));
        auto A22_right = A22(_,_(N/2+1,   N));

        const auto A12_left  = A12(_,_(    1, N/2));
        const auto A12_right = A12(_,_(N/2+1,   N));

        auto gemm0 = std::thread([&]{ A22_left  -= A21*A12_left;  });
        auto gemm1 = std::thread([&]{ A22_right -= A21*A12_right; });
        gemm0.join();
        gemm1.join();

        /*
        const IndexType M = A22.numRows();

        auto A22_top = A22(_(    1, M/2), _);
        auto A22_bot = A22(_(M/2+1,   M), _);

        const auto A21_top = A21(_(    1, M/2), _);
        const auto A21_bot = A21(_(M/2+1,   M), _);

        auto gemm0 = std::thread([&]{ A22_top -= A21_top*A12; });
        auto gemm1 = std::thread([&]{ A22_bot -= A21_bot*A12; });
        gemm0.join();
        gemm1.join();
        */
    }
    return 0;
}

} // namespace flens

#endif // LU_LU_UNBLK_MT_H
