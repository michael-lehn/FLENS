#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    typedef complex<double>             Complex;
    const Complex                       I(0,1);

    GeMatrix<FullStorage<Complex> >     A(4, 4);
    DenseVector<Array<Complex> >        b(4);
    DenseVector<Array<int> >            piv(4);

///
/// Setup matrix $A$.  Note that `2+I` or `2*I` would not work.  That's because
/// the STL does not define the operation `int+complex<double>` or
/// `int*complex<double>`.
///
    A =   2,    3,  -I,    0,
         -6,   -5,   0, 2.+I,
       2.*I,   -5,   6,   -6,
          4, 6.*I,   2,   -3;

    b = 20,
       -33.*I,
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
