#include <iostream>
///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx
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
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    ///
    ///  Define an underscore operator for convenient matrix slicing
    ///
    const Underscore<IndexType> _;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType m = 4,
                    n = 5;

    Matrix            Ab(m, n);
    IndexVector       piv(m);

    Ab =  2,   3,  -1,   0,  20,
         -6,  -5,   0,   2, -33,
          2,  -5,   6,  -6, -43,
          4,   6,   2,  -3,  49;

    cout << "Ab = " << Ab << endl;

    ///
    /// Compute the $LU$ factorization with __lapack::trf__
    ///
    lapack::trf(Ab, piv);

    ///
    ///  Solve the system of linear equation $Ax =B$ using __blas::sm__
    ///
    const auto A = Ab(_,_(1,m));
    auto       B = Ab(_,_(m+1,n));

    blas::sm(Left, NoTrans, T(1), A.upper(), B);

    cout << "X = " << B << endl;
}

///
///  :links: __lapack::trf__ -> file:flens/lapack/ge/trf.h
///          __blas::sm__ -> file:flens/blas/level3/sm.h
///
