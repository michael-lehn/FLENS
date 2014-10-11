#include <cxxstd/iostream.h>
#include <flens/examples/my_symatrix.h>
#include <flens/examples/cg.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int   n = 5;

    MySyMatrix  A;

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
