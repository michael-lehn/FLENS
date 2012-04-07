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
    ///
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef DenseVector<Array<T> >      DenseVector;

    ///
    ///  Define an underscore operator for convenient matrix slicing
    ///
    typedef GeMatrix::IndexType   IndexType;
    const Underscore<IndexType>   _;

    const IndexType m = 4,
                    n = 5;

    GeMatrix  Ab(m, n);

    Ab =  2,   3,  -1,   0,   20,
         -6,  -5,   0,   2,  -33,
          2,  -5,   6,  -6,  -43,
          4,   6,   2,  -3,   49;

    cout << "Ab = " << Ab << endl;

    ///
    ///  `tau` will contain the scalar factors of the elementary reflectors
    ///  (see __dgeqrf__ for details).  Vector `work` is used as workspace
    ///  and, if empty, gets resized by __lapack::qrf__ to optimal size.
    ///
    DenseVector tau(min(m,n));
    DenseVector work;

    ///
    ///  Compute the $QR$ factorization of $A$ with __lapack::qrf__ the
    ///  FLENS port of LAPACK's __dgeqrf__.
    ///
    lapack::qrf(Ab, tau, work);
    cout << "Ab = " << Ab << endl;
    cout << "tau = " << tau << endl;

    ///
    ///  Solve the system of linear equation $Ax =B$ using the triangular
    ///  solver __blas::sm__ which is the FLENS implementation of the BLAS
    ///  routine `trsm`.  Note that `A.upper()` creates a triangular view.
    ///
    const auto A = Ab(_,_(1,m));
    auto       B = Ab(_,_(m+1,n));

    blas::sm(Left, NoTrans, 1, A.upper(), B);

    cout << "X = " << B << endl;
}

///
///  :links: __lapack::qrf__ -> file:flens/lapack/impl/qrf.h
///          __dgeqrf__      -> file:cxxlapack/netlib/lapack/dgeqrf.f
///          __blas::sm__    -> file:flens/blas/level3/sm.h
///
