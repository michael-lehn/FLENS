#ifndef LU_LU_PANEL_MT_H
#define LU_LU_PANEL_MT_H 1

#ifndef LU_PANEL_BS
#define LU_PANEL_BS 32
#endif

#ifndef LU_PANEL_M_BS
#define LU_PANEL_M_BS 128
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/threadpool.h>
#include <flens/examples/lu/lu_unblk_with_operators.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_panel(ThreadPool &tp, GeMatrix<MA> &A, DenseVector<VP> &p)
{
    typedef typename GeMatrix<MA>::ElementType  T;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const IndexType  m = A.numRows();
    const IndexType  n = A.numCols();
    const IndexType  mn = std::min(m, n);

    const T One(1);

    const Underscore<IndexType>  _;

    IndexType info = 0;

    for (IndexType i=1; i<=mn; i+=LU_PANEL_BS) {
        const IndexType bs = std::min(std::min(LU_PANEL_BS, m-i+1), n-i+1);

        // Partitions of A
        auto A_0 = A(_(   i,      m), _(   1,    i-1));
        auto A_1 = A(_(   i,      m), _(   i, i+bs-1));

        auto A11 = A(_(   i, i+bs-1), _(   i, i+bs-1));
        auto A21 = A(_(i+bs,      m), _(   i, i+bs-1));
        auto A12 = A(_(   i, i+bs-1), _(i+bs,      n));
        auto A22 = A(_(i+bs,      m), _(i+bs,      n));

        // Part of the pivot vector for rows of A11
        auto p_  = p(_(i,m));

        // Compute LU factorization of A11
        info = lu_unblk(A_1, p_);

        if (info) {
            // All values in column info of A11 are *exactly* zero.
            return info+i-1;
        }

        // Apply permutation to A_0
        auto task = [=]() mutable
                    {
                        apply_perm_inv(p(_(i,i+bs-1)), A_0);
                    };
        tp.add(task);


        for (IndexType l=i+bs; l<=n; l+= LU_PANEL_BS) {
            const IndexType nb  = std::min(LU_PANEL_BS, n-l+1);
            auto            A_L = A(_(   i,      m), _(   l, l+nb-1));
            auto            A1L = A(_(   i, i+bs-1), _(   l, l+nb-1));

            // Apply permutation to A_L
            // Use triangular solver for A1L = A11.lowerUnit()*A1L
            auto task = [=]() mutable
                        {
                            apply_perm_inv(p(_(i,i+bs-1)), A_L);
                            blas::sm(Left, NoTrans, One, A11.lowerUnit(), A1L);
                        };
            tp.add(task);
        }
        tp.join();

        // Update p
        p_ += i-1;

        for (IndexType l=i+bs; l<=n; l+= LU_PANEL_BS) {
            const IndexType nb  = std::min(LU_PANEL_BS, n-l+1);
            auto            A_L = A(_(   i,      m), _(   l, l+nb-1));
            auto            A1L = A(_(   i, i+bs-1), _(   l, l+nb-1));

            for (IndexType k=i+bs; k<=m; k+=LU_PANEL_M_BS) {
                const IndexType mb  = std::min(LU_PANEL_M_BS, m-k+1);
                const auto      AK1 = A(_(k,k+mb-1), _(   i, i+bs-1));
                auto            AKL = A(_(k,k+mb-1), _(   l, l+nb-1));


                // Update AKL with matrix-product AKL = AKL - AK1*A1L
                auto task = [=]() mutable
                            {
                                blas::mm(NoTrans, NoTrans,
                                         -One, AK1, A1L,
                                         One, AKL);
                            };
                tp.add(task);
            }
        }
        tp.join();
    }
    return 0;
}

} // namespace flens

#endif // LU_LU_PANE_MTL_H
