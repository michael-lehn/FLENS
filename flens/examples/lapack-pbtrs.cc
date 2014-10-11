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
    ///  We initialize the positiv definite $n \times n$ tridiagonal matrix.
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
    /// We factorize $A$ with __lapack::pbtrf__
    ///
    /// :links: __lapack::pbtrf__ -> doc:flens/lapack/pb/pbtrf
    lapack::pbtrf(A);

    ///
    /// In the upper part of $A$ now the triangular factor $L^T$ is stored.
    ///
    cout << "L^T = " << A.general().upper() << endl;

    ///
    /// Solve it with __lapack::pbtrs__
    ///
    /// :links: __lapack::pbtrs__ -> doc:flens/lapack/pb/pbtrs
    lapack::pbtrs(A, b);

    cout << "x = " << b << endl;
}
