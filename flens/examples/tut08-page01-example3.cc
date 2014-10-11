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
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;

///
/// We want a matrix with kl sub- and ku superdiagonals,
/// thus we have to declare it with kl sub- and kl+ku superdiagonals.
///

    Matrix         Ab(n, n, 2, 1+2);
    Vector         b(n), x(n);
    IndexVector    pivots(n);

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
    /// We compute $x = A*b$.
    ///
    x = Ab*b;

    cout << "x = " << x << endl;

    ///
    /// And solve  $A*x = b$.
    ///
    lapack::sv(Ab, pivots, x);

    cout << "b = " << x << endl;

    return 0;
}
