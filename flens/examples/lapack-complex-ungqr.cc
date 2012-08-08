#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    typedef complex<double>             Complex;
    const Complex                       I(0,1);


    GeMatrix<FullStorage<Complex> >     A(4,4), Q;
    DenseVector<Array<Complex> >        b(4);
    DenseVector<Array<Complex> >        tau;
    //DenseVector<Array<Complex> >      work;

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;
    A *= I;

    b = 20,
       -33,
       -43,
        49;
    b *= I;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the factorization $A = QR$.  Note that the workspace gets created
/// implicitly and temporarily.  So you might not want to do this inside a loop.
///
    lapack::qrf(A, tau);
    //lapack::qrf(A, tau, work);

///
/// Explicitly setup $Q$.
///
    Q = A;
    lapack::ungqr(Q, tau);
    //lapack::orgqr(Q, tau, work);

    cout << "Q = " << Q << endl;

///
/// Compute $\tilde{b} = Q^T b$.
///
    DenseVector<Array<Complex> >  x;
    x = conjTrans(Q)*b;

///
/// Solve $R x = \tilde{b}$.  Vector $b$ gets overwritten with $x$.
///
    const auto R = A.upper();
    blas::sv(NoTrans, R, x);

    cout << "x = " << x << endl;
}
