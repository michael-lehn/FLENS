#include <flens/flens.cxx>
#include <iostream>

#ifdef WITH_OPERATORS
#   include <flens/examples/lu/lu_unblk_with_operators.h>
#else
#   include <flens/examples/lu/lu_unblk_with_flensblas.h>
#endif

#include <flens/examples/lu/apply_perm.h>
#include <flens/examples/lu/timer.h>

using namespace flens;
using namespace std;

int
main()
{
    const int m = 2000;
    const int n = 2000;

    GeMatrix<FullStorage<double> >  A(m,n), A_org;
    DenseVector<Array<int> >        p(m);
    Underscore<int>                 _;

    fillRandom(A);

    // store original A for checking the factorization later
    A_org = A;

    /// Uncomment for playing with small dimensions:
    //cout << "A = " << A << endl;

    double t0 = ATL_walltime();
    if (lu_unblk(A, p) != 0) {
        cout << "Matrix is singular or close to singular" << endl;
        return 1;
    }
    t0 = ATL_walltime() - t0;

    cout << "Time elapsed: " << t0 << endl;

    /// and
    //cout << "L = " << A.lowerUnit() << endl;
    //cout << "U = " << A.upper() << endl;
    //cout << "p = " << p << endl;

    cout << "Checking results:" << endl;
    GeMatrix<FullStorage<double> >  LU(m,n);

    // blas::copy(NoTrans, A.upper(), LU.upper());
    LU.upper() = A.upper();

    // blas::mm(Left, NoTrans, 1.0, A(_(1,4),_(1,4)).lowerUnit(), LU);
    LU = A(_(1,m),_(1,m)).lowerUnit() * LU;

    // LU = P*LU
    apply_perm(p, LU);

    // blas::axpy(NoTrans, -1.0, LU, A_org);
    A_org -= LU;

    cout << "|| A-P*L*U ||_1 = " << blas::asum(A_org) << endl;
    return 0;
}
