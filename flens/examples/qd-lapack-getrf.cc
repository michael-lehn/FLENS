#include <iostream>

#include <qd/qd_real.h>
#include <qd/fpu.h>

///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef qd_real   T;

int
main()
{
    unsigned int old_cw;
    fpu_fix_start(&old_cw);


    T a = 2;
    T b = 1;
    T c = 3;

    cout << "min(a,b,c) = " << flens::min(a,b,c) << endl;
    cout << "max(a,b,c) = " << flens::max(a,b,c) << endl;
    cout << a[0] << endl;
    cout << a.x[0] << endl;

    cout << "sizeof(T) = " << sizeof(T) << endl;
    cout << "numeric_limits<double>::epsilon() = "
         << numeric_limits<double>::epsilon()
         << endl;
    cout << "numeric_limits<T>::epsilon() = "
         << numeric_limits<T>::epsilon()
         << endl;
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

    blas::sm(Left, NoTrans, 1, A.upper(), B);

    cout << "X = " << B << endl;

    fpu_fix_end(&old_cw);
}

///
///  :links: __lapack::trf__ -> file:flens/lapack/gesv/trf.h
///          __blas::sm__ -> file:flens/blas/level3/sm.h
///
