#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    ///
    ///  Define some convenient typedefs for the matrix/vector types
    ///  of our system of linear equations.
    ///
    typedef SbMatrix<BandStorage<T> >           SymmetricBandMatrix;
    typedef DenseVector<Array<T> >              Vector;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef SymmetricBandMatrix::IndexType      IndexType;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType n = 4;

    ///
    ///  We allocate a symmetric band matrix with one off-diagonal.
    ///
    SymmetricBandMatrix     A(n, Upper, 1);

    Vector                  b(n);

    ///
    /// We initialize the positiv definite $n \times n$ tridiagonal matrix.
    /// Note that currently the method `diag(d)` is only provided for matrices
    /// of type `GbMatrix` not `SbMatrix`.  So we access the method via a
    /// `GbMatrix` view.  If frequently needed we can provide the method also
    /// through the `SbMatrix` interface directly.
    ///
    A.general().diag( 1) = -1;
    A.general().diag( 0) =  2;

    b =  1,
         2,
         3,
         4;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

    ///
    /// And solve it with __lapack::pbsv__
    ///
    /// :links: __lapack::pbsv__ -> doc:flens/lapack/pb/pbsv
    lapack::pbsv(A, b);

    cout << "(L\\R) = " << A << endl;

    cout << "x = " << b << endl;
}
