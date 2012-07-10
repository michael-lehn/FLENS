#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4, 4);
    DenseVector<Array<double> >      b(4);
    DenseVector<Array<int> >         piv(4);

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

    b = 20,
       -33,
       -43,
        49;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the $PLU = A$ factorization with __lapack::trf__
///
    lapack::trf(A, piv);

///
/// Solve the system of linear equation $PLUx = b$ using __lapack::trs__.
/// Vector $b$ gets overwritten with vector $x$.
///
    lapack::trs(NoTrans, A, piv, b);
    cout << "X = " << b << endl;
}

///
///  :links: __lapack::trf__ -> file:flens/lapack/impl/trf.h
///          __lapack::trs__ -> file:flens/lapack/impl/trs.h
///
