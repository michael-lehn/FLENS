#ifndef LU_LU_PANEL_BLK_H
#define LU_LU_PANEL_BLK_H 1

#ifndef LU_BS
#define LU_BS 256
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_panel.h>
#include <flens/examples/lu/timer.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_panel_blk(GeMatrix<MA> &A, DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = std::min(m, n);

    const T One(1);

    const Underscore<IndexType>  _;

    IndexType info = 0;

    //double T1 = 0;
    //double T2 = 0;

    for (IndexType i=1; i<=mn; i+=LU_BS) {
        const IndexType bs = std::min(std::min(LU_BS, m-i+1), n-i+1);

        // Partitions of A
        auto A_0 = A(_(   i,      m), _(   1,    i-1));
        auto A_1 = A(_(   i,      m), _(   i, i+bs-1));
        auto A_2 = A(_(   i,      m), _(i+bs,      n));

        auto A11 = A(_(   i, i+bs-1), _(   i, i+bs-1));
        auto A21 = A(_(i+bs,      m), _(   i, i+bs-1));
        auto A12 = A(_(   i, i+bs-1), _(i+bs,      n));
        auto A22 = A(_(i+bs,      m), _(i+bs,      n));

        // Part of the pivot vector for rows of A11
        auto p_  = p(_(i,m));

        //double t0 = ATL_walltime();

        // Compute LU factorization of A11
        info = lu_panel(A_1, p_);


        if (info) {
            // All values in column info of A11 are *exactly* zero.
            return info+i-1;
        }

        // Apply permutation to A10 and A12
        apply_perm_inv(p(_(i,i+bs-1)), A_0);
        apply_perm_inv(p(_(i,i+bs-1)), A_2);
        //double t1 = ATL_walltime();

        // Update p
        p_ += i-1;

        // Use triangular solver for A12 = A11.lowerUnit()*A12
        blas::sm(Left, NoTrans, One, A11.lowerUnit(), A12);

        // Update A22 with matrix-product A22 = A22 - A21*A12
        blas::mm(NoTrans, NoTrans, -One, A21, A12, One, A22);

        //double t2 = ATL_walltime();

        //T1 += t1 - t0;
        //T2 += t2 - t1;
    }
    //std::cout << "T1 = " << T1 << ", T2 = " << T2 << std::endl;
    return 0;
}

} // namespace flens

#endif // LU_LU_PANEL_BLK_H
