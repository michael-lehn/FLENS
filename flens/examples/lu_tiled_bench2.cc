#include <flens/flens.cxx>
#include <iostream>

#include <flens/examples/lu/apply_perm.h>
#include <flens/examples/lu/check_lu.h>
#include <flens/examples/lu/lu_tiled_mt2.h>
#include <flens/examples/lu/scheduler2.h>
#include <flens/examples/lu/timer.h>

#ifndef NUM_THREADS
#define NUM_THREADS 2
#endif

using namespace flens;
using namespace std;

typedef GeMatrix<FullStorage<double> >  DGeMatrix;
typedef TiledCopy<DGeMatrix::View>      DTiledMatrix;
typedef DenseVector<Array<double> >     DDenseVector;
typedef DenseVector<Array<int> >        IDenseVector;


int
main()
{
    const int m0    = 256;
    const int n0    = 256;

    const int mInc  = 128;
    const int nInc  = 128;

    const int m1    = 6000;
    const int n1    = 6000;
    const int mn1   = std::min(m1, n1);

    const int numIt = 10;

    Scheduler scheduler(NUM_THREADS);
    int bs = 256;

    DGeMatrix           A_buffer(m1,n1), A_org_buffer(m1,n1);
    IDenseVector        p_buffer(mn1);

    Underscore<int>     _;

    for (int m=m0, n=n0; m<=m1 && n<=n1; m+=mInc, n+=nInc) {

        auto A     = A_buffer(_(1,m),_(1,n));
        auto A_org = A_org_buffer(_(1,m),_(1,n));
        auto p     = p_buffer(_(1,std::min(m,n)));

        if (std::min(m,n) < 256) {
            bs = 64;
        } else if (std::min(m,n) < 512) {
            bs = 128;
        } else if (std::min(m,n) < 3*1024) {
            bs = 256;
        } else {
            bs = 512;
        }

        double t_elapsed = 0.0;
        double err_nrm   = 0.0;


        for (int it=0; it<numIt; ++it) {
            fillRandom(A);
            A_org = A;

            DTiledMatrix        A_(A, bs);

            double t0 = ATL_walltime();
            int info  = lu_tiled(scheduler, A_, p);
            t0 = ATL_walltime() - t0;

            A_.untile(A);

            if (info!=0) {
                --it;
            } else {
                //err_nrm   += check_LU(A_org, A, p);
                t_elapsed += t0;
            }
        }
        err_nrm   /= numIt;
        t_elapsed /= numIt;

        std::string report = format("%6d %6d %9.3f %10.2e",
                                    m, n, t_elapsed, err_nrm);
        cout << report << endl;
    }
}
