#include <cmath>
#include <iostream>
#include <flens/flens.cxx>

#ifndef USE_CG_BLAS
#include <flens/examples/cg.h>
#else
#include <flens/examples/cg_blas.h>
#endif

using namespace flens;
using namespace std;

int
main()
{
    const int n = 5;

///
/// Let the diagonal of A have values of `2` and the off-diagonal `-1`.
///
    SyMatrix<FullStorage<double> >   A(n, Upper);
    A.diag(0) = 2;
    A.diag(1) = -1;

    cout << "A = " << A << endl;
///
/// Vector `z` is the control solution.
///
    DenseVector<Array<double> > x(n), z(n), b(n);
    z = 1;
    b = A*z;

///
/// Solve $Ax=b$
///
    cg(A, b, x);
    cout << "x = " << x << endl;

///
/// Get the infinity norm of the error
///
    x -= z;
    cout << "max abs error = " << abs(x(blas::iamax(x))) << std::endl;
}
