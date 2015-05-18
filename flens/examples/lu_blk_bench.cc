#include <flens/flens.cxx>
#include <algorithm>
#include <iostream>

#ifdef WITH_OPERATORS
#   include <flens/examples/lu/lu_blk_with_operators.h>
#else
#   include <flens/examples/lu/lu_blk_with_flensblas.h>
#endif

#include <flens/examples/lu/check_lu.h>
#include <flens/examples/lu/timer.h>
#include <flens/examples/lu/format.h>

#define POW3(x) (x)*(x)*(x)

using namespace flens;
using namespace std;

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

    GeMatrix<FullStorage<double> >  A_buffer(m1,n1), A_org_buffer(m1,n1);
    DenseVector<Array<int> >        p_buffer(mn1);

    Underscore<int>                 _;


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
            int info  = lu_blk(A, p);
            t0 = ATL_walltime() - t0;

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
    return 0;
}
