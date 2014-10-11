#include <cxxstd/cmath.h>
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>
#include <flens/examples/cg.h>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 5;

///
/// Let the diagonal of A have values of `2` and the off-diagonal `-1`.
///
    SyCoordMatrix<CoordStorage<double> >    A(n, Upper);
    for (int i=1; i<=n; ++i) {
        A(i,i) += 2;
        if (i<n) {
            A(i,i+1) -= 1;
        }
    }

///
/// Compress A and stor it in B
///
    SyCRSMatrix<CRS<double> >   B = A;

///
/// Vector `z` is the control solution.
///
    DenseVector<Array<double> > x(n), z(n), b(n);
    z = 1;
    b = B*z;

///
/// Solve $Bx=b$
///
    cg(B, b, x);
    cout << "x = " << x << endl;

///
/// Get the infinity norm of the error
///
    x -= z;
    cout << "max abs error = " << abs(x(blas::iamax(x))) << std::endl;
}
