#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef complex<double>   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const Underscore<IndexType> _;

    ///
    ///  Set up the baby problem ...
    ///
    const IndexType m = 4,
                    n = 5;

    Matrix            Ab(m, n);
    IndexVector       piv(m);

    Ab = T(1, 0), T( 1,-1), T(  2,20), T( 1, 0),
         T(0, 1), T( 1, 2), T(-10, 8), T(-1, 1),
         T(0,-1), T(-1, 1), T( 40, 6), T(-2,-1);

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
    cout << "piv = " << piv << endl;
}

///
///  :links: __lapack::trf__ -> file:flens/lapack/ge/trf.h
///          __blas::sm__ -> file:flens/blas/level3/sm.h
///
