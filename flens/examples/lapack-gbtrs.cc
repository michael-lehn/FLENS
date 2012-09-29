#include <iostream>
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
    typedef GbMatrix<BandStorage<T> >           BandMatrix;
    typedef DenseVector<Array<T> >              Vector;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef BandMatrix::IndexType               IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType n = 4;

    ///
    ///  We allocate a band matrix with one sub-diagonal and two superdiagonals.
    ///  Note, if you compute the LU factorization of a tridiagonal matrix then
    ///  an extra superdiagonal is needed!
    ///
    BandMatrix     A(n,n,1,2);

    Vector         b(n);
    IndexVector    piv(n);

    ///
    ///  We initialize the $n \times n$ tridiagonal matrix.
    ///
    A.diag( 1) = -1;
    A.diag( 0) =  2;
    A.diag(-1) = -3;

    b = 4,
       -7,
        0,
        1;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

    ///
    /// Compute the $LU$ factorization with __lapack::trf__
    ///
    /// :links: __lapack::trf__ -> doc:flens/lapack/gb/trf
    lapack::trf(A, piv);

    cout << "(L\\R) = " << A << endl;

    ///
    /// Solve the factorized problem with __lapack::trs__
    ///
    /// :links: __lapack::trs__ -> doc:flens/lapack/gb/trs
    lapack::trs(NoTrans, A, piv, b);

    cout << "x = " << b << endl;
}

