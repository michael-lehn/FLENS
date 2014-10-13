#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GbMatrix<BandStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;

    const IndexType n = 4;

///
/// The first two indices define the dimensions, the two last ones the
/// the number of sub- and superdiagonals.
///

    Matrix         Ab(n, n, 2, 1);
    Vector         b(n), x(n);

    Ab.diag(-2) =  2, 3;
    Ab.diag(-1) = -6, 5, -4;
    Ab.diag( 0) =  4, 3,  5, -1;
    Ab.diag( 1) =  1, 2,  7;

    b = 20,
       -33,
       -43,
        49;

    cout << "Ab = " << Ab << endl;
    cout << "b  = " << b << endl;

    ///
    /// The function lower() creates a triangular banded matrix.
    /// We compute $x = A_{lower}^T b$.
    ///
    x = transpose(Ab.lower())*b;

    cout << "x = " << x << endl;

    ///
    /// And solve  $x = A_{lower}^T b$.
    ///
    blas::sv(Trans, Ab.lower(), x);

    cout << "b = " << x << endl;

    return 0;
}
