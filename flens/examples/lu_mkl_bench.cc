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
typedef DenseVector<Array<double> >     DDenseVector;
typedef DenseVector<Array<long> >        IDenseVector;


extern "C" {

#define MKL_INT long

void
dgetrf(const MKL_INT* m, const MKL_INT* n, double* a, const MKL_INT* lda,
       MKL_INT* ipiv, MKL_INT* info );

}

int
main()
{
    const int m0    = 256;
    const int n0    = 256;

    const int mInc  = 256;
    const int nInc  = 256;

    const int m1    = 10000;
    const int n1    = 10000;
    const int mn1   = std::min(m1, n1);

    const int numIt = 10;

    DGeMatrix           A_buffer(m1,n1), A_org_buffer(m1,n1);
    IDenseVector        p_buffer(mn1);

    Underscore<int>     _;

    for (int m=m0, n=n0; m<=m1 && n<=n1; m+=mInc, n+=nInc) {

        auto A     = A_buffer(_(1,m),_(1,n));
        auto A_org = A_org_buffer(_(1,m),_(1,n));
        auto p     = p_buffer(_(1,std::min(m,n)));

        double t_elapsed = 0.0;
        double err_nrm   = 0.0;


        for (int it=0; it<numIt; ++it) {
            fillRandom(A);
            A_org = A;


            double t0 = ATL_walltime();
            long info;
            long M   = m;
            long N   = m;
            long ldA = A.leadingDimension();
            dgetrf(&M, &N, A.data(), &ldA, p.data(), &info);
            t0 = ATL_walltime() - t0;

            if (info!=0) {
                --it;
            } else {
                err_nrm   += check_LU(A_org, A, p);
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
