#ifndef LU_LU_TILED_MT_H
#define LU_LU_TILED_MT_H 1

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cstdio>
#include <vector>
#include <flens/flens.h>

#include <flens/examples/lu/apply_perm_inv.h>
#include <flens/examples/lu/lu_blk_with_operators.h>
#include <flens/examples/lu/lu_tiled_keys2.h>
#include <flens/examples/lu/scheduler2.h>
#include <flens/examples/lu/tile.h>

namespace flens {

template <typename MA, typename VP>
typename DenseVector<VP>::IndexType
lu_tiled(Scheduler &scheduler, TiledCopy<MA> &A, DenseVector<VP> &p)
{
    typedef typename TiledCopy<MA>::ElementType  T;
    typedef typename TiledCopy<MA>::IndexType    IndexType;
    typedef Scheduler::Task                      Task;

    const IndexType  m  = A.numTileRows();
    const IndexType  n  = A.numTileCols();
    const IndexType  bs = A.blockSize();
    const IndexType  mn = std::min(m, n);

    const T One(1);
    const Underscore<IndexType>  _;

    Task task;

    for (IndexType i=1; i<=mn; ++i) {
        auto A_ii = A(i,i);
        auto p_   = p(_(bs*(i-1)+1,bs*(i-1)+A_ii.numRows()));

        //  Factorize A(i,i)
        task = [=, &scheduler]() mutable
               {
                   auto info = lu_blk(A_ii, p_);
                   if (info!=0) {
                       scheduler.abort(info+bs*(i-1));
                   }
               };
        if (i>1) {
            scheduler.addArc(keyUpdateA(mn,m,n,i,i,i-1), keyLU(i));
        }
        scheduler.add(keyLU(i), task);

        //  Apply permutation to i-th tile row
        for (IndexType j=1; j<=n; ++j) {
            if (i!=j) {
                auto A_ij = A(i,j);
                task = [=]() mutable
                       {
                           apply_perm_inv(p_, A_ij);
                       };
                scheduler.add(keyPtA(mn,m,i,j), task, { keyLU(i) });

                // dependency for task: "Update Permutation Vector p"
                scheduler.addArc(keyPtA(mn, m, i, j), keyUpdateP(mn,m,n,i));
            }
        }

        //  Update Permutation Vector p
        task = [=]() mutable
               {
                   p_ += bs*(i-1);
               };
        if (i>1) {
            scheduler.addArc(keyUpdateP(mn,m,n,i-1), keyUpdateP(mn,m,n,i));
        }
        scheduler.add(keyUpdateP(mn,m,n,i), task);

        //  Apply U_ii^{-1} to i-th tile column
        for (IndexType k=i+1; k<=m; ++k) {
            auto A_ki = A(k,i);
            task = [=]() mutable
                   {
                       auto U_ii = A_ii.upper();
                       blas::sm(Right, NoTrans, One, U_ii, A_ki);
                   };
            if (i>1) {
                scheduler.addArc(keyUpdateA(mn,m,n,k,i,i-1),
                                 keyApplyU(mn,m,n,k,i));
            }
            scheduler.add(keyApplyU(mn,m,n,k,i), task, { keyLU(i) });
        }

        //  Apply L_ii^{-1} to i-th tile row
        for (IndexType l=i+1; l<=n; ++l) {
            auto A_il = A(i,l);
            task = [=]() mutable
                   {
                       auto L_ii = A_ii.lowerUnit();
                       blas::sm(Left, NoTrans, One, L_ii, A_il);
                   };
            if (i>1) {
                scheduler.addArc(keyUpdateA(mn,m,n,i,l,i-1),
                                 keyApplyL(mn,m,n,i,l));
            }
            scheduler.add(keyApplyL(mn,m,n,i,l), task, { keyPtA(mn,m,i,l) });
        }

        //  Update remaining blocks
        for (IndexType k=i+1; k<=m; ++k) {
            const auto A_ki = A(k,i);
            for (IndexType l=i+1; l<=n; ++l) {
                const auto A_il = A(i,l);
                auto       A_kl = A(k,l);

                task = [=]() mutable
                       {
                           A_kl -= A_ki * A_il;
                       };
                scheduler.add(keyUpdateA(mn,m,n,k,l,i), task,
                              { keyApplyU(mn,m,n,k,i), keyApplyL(mn,m,n,i,l) });
            }
        }
    }
    return scheduler.join();
}

} // namespace flens

#endif // LU_LU_TILED_MT_H
