#ifndef LU_LU_PANEL_BLK_MT_H
#define LU_LU_PANEL_BLK_MT_H 1

#ifndef LU_BS
#define LU_BS 256
#endif

#ifndef LU_N_BS
#define LU_N_BS 4096
#endif

#ifndef LU_M_BS
#define LU_M_BS 386
#endif

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_panel_mt.h>
#include <flens/examples/lu/threadpool.h>
#include <flens/examples/lu/timer.h>
#include <flens/examples/lu/lu_unblk_with_operators.h>

namespace flens {

template <typename MA, typename VP>
typename GeMatrix<MA>::IndexType
lu_panel_blk(ThreadPool &tp, GeMatrix<MA> &A, DenseVector<VP> &p)
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
    //double T3 = 0;

    for (IndexType i=1; i<=mn; i+=LU_BS) {
        const IndexType bs = std::min(std::min(LU_BS, m-i+1), n-i+1);

        // Panel of A
        auto A_0 = A(_(   i,      m), _(   1,    i-1));
        auto A_1 = A(_(   i,      m), _(   i, i+bs-1));
        //auto A_2 = A(_(   i,      m), _(i+bs,      n));

        // Part of the pivot vector for rows of A11
        auto p_  = p(_(i,m));

        //double t0 = ATL_walltime();

        // Compute LU factorization of A11
        info = lu_panel(tp, A_1, p_);
        //info = lu_panel(A_1, p_);

        //double t1 = ATL_walltime();

        if (info) {
            // All values in column info of A11 are *exactly* zero.
            return info+i-1;
        }

        // Apply permutation to A10 and A12
        //apply_perm_inv(p(_(i,i+bs-1)), A_0);
        //apply_perm_inv(p(_(i,i+bs-1)), A_2);

        for (IndexType l=1; l<i; l+=LU_BS) {
            const IndexType nb  = std::min(LU_BS, i-l);
            auto            A_L = A(_(   i,      m), _(   l, l+nb-1));

            auto task = [=]() mutable
                        {
                            apply_perm_inv(p(_(i,i+bs-1)), A_L);
                        };
            tp.add(task);
        }
        tp.join();

        //double t2 = ATL_walltime();

        const auto A11 = A(_(   i, i+bs-1), _(   i, i+bs-1));
        for (IndexType l=i+bs; l<=n; l+=LU_BS/2) {
            const IndexType nb  = std::min(LU_BS/2, n-l+1);
            auto            A_L = A(_(   i,      m), _(   l, l+nb-1));
            auto            A1L = A(_(   i, i+bs-1), _(   l, l+nb-1));
            auto            A2L = A(_(i+bs,      m), _(   l, l+nb-1));
            auto            A21 = A(_(i+bs,      m), _(   i, i+bs-1));

            // Use triangular solver for A1L = A11.lowerUnit()*A1L
            auto task = [=]() mutable
                        {
                            apply_perm_inv(p(_(i,i+bs-1)), A_L);
                            blas::sm(Left, NoTrans, One, A11.lowerUnit(), A1L);
                            //blas::mm(NoTrans, NoTrans, -One, A21, A1L, One, A2L);
                        };
            tp.add(task);
        }
        tp.join();

        // Update p
        p_ += i-1;


        for (IndexType l=i+bs; l<=n; l+=LU_N_BS) {
            const IndexType nb  = std::min(LU_N_BS, n-l+1);
            const auto      A1L = A(_(   i, i+bs-1), _(   l, l+nb-1));

            for (IndexType k=i+bs; k<=m; k+=LU_M_BS) {
                const IndexType mb  = std::min(LU_M_BS, m-k+1);
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

        //double t3 = ATL_walltime();
        //T1 += t1 - t0;
        //T2 += t2 - t1;
        //T3 += t3 - t2;

    }
    //std::cout << "T1 = " << T1
    //          << ", T2 = " << T2
    //          << ", T3 = " << T3 << std::endl;
    return 0;
}

} // namespace flens

#endif // LU_LU_PANEL_BLK_MT_H
