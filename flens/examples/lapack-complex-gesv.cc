#include <complex>
#include <iostream>
///
/// With header __flens.cxx__ all of FLENS gets included.
///
/// :links:  __flens.cxx__ -> file:flens/flens.cxx
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

///
/// Elements are now defined to be *complex*
///
typedef complex<double>   T;

int
main()
{
    ///
    /// Like before we define convenient matrix/vector types ...
    ///
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    ///
    /// Then we setup another toy problem
    ///
    const IndexType n = 3;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    A = T(1, 0), T( 1,-1), T(  2,20),
        T(0, 1), T( 1, 2), T(-10, 8),
        T(0,-1), T(-1, 1), T( 40, 6);

    b = T( 1, 0),
        T(-1, 1),
        T(-2,-1);

    cerr << "A = " << A << endl;
    cerr << "b = " << b << endl;

    ///
    /// And solve it with __lapack::sv__
    ///
    /// :links: __lapack::sv__ -> doc:flens/lapack/ge/sv
    lapack::sv(A, piv, b);

    cerr << "A = " << A << endl;
    cerr << "piv = " << piv << endl;
    cerr << "x = " << b << endl;
}

