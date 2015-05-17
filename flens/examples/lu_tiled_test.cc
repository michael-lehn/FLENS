#include <flens/flens.cxx>
#include <iostream>

#include <flens/examples/lu/apply_perm.h>
#include <flens/examples/lu/check_lu.h>
#include <flens/examples/lu/lu_tiled_mt.h>
#include <flens/examples/lu/scheduler.h>
#include <flens/examples/lu/timer.h>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >  DGeMatrix;
    typedef TiledCopy<DGeMatrix>            DTiledMatrix;
    typedef DenseVector<Array<double> >     DDenseVector;
    typedef DenseVector<Array<int> >        IDenseVector;

    Scheduler scheduler(2);

    const int bs = 256;
    const int m  = 4000;
    const int n  = 4000;

    DGeMatrix           A(m,n), A_org;
    IDenseVector        p(m);
    Underscore<int>     _;

    // Init A and store a copy for checking the factorization later
    fillRandom(A);
    A_org = A;

    double t0, t1;

    // Create a new matrix and initialize it with the tiled format
    t0 = ATL_walltime();

    DTiledMatrix        A_(A, bs);

    t1 = ATL_walltime();

    // Call the tile algorithm
    if (lu_tiled(scheduler, A_, p) != 0) {
        cout << "Matrix is singular or close to singular" << endl;
        return 1;
    }
    t1 = ATL_walltime()-t1;

    // Copy back the tiled matrix
    A_.untile(A);

    t0 = ATL_walltime()-t0;

    cout << "total time elapsed:          " << t0 << endl;
    cout << "time elapsed for 'lu_tiled': " << t1 << endl;
    cout << "time elapsed for conversion: " << t0-t1 << endl;

    // Check result
    cout << "Checking results:" << endl;
    cout << "|| A-P*L*U ||_1 = " << check_LU(A_org, A, p) << endl;
    return 0;
}
