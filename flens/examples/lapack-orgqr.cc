#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    GeMatrix<FullStorage<double> >     A(4,4), Q;
    DenseVector<Array<double> >        b(4);
    DenseVector<Array<double> >        tau;
    //DenseVector<Array<double> >      work;

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
/// Compute the factorization $A = QR$.  Note that the workspace gets created
/// implicitly and temporarily.  So you might not want to do this inside a loop.
///
    lapack::qrf(A, tau);
    //lapack::qrf(A, tau, work);

///
/// Explicitly setup $Q$.
///
    Q = A;
    lapack::orgqr(Q, tau);
    //lapack::orgqr(Q, tau, work);

    cout << "Q = " << Q << endl;

///
/// Compute $\tilde{b} = Q^T b$.
///
    DenseVector<Array<double> >  _b;
    _b = transpose(Q)*b;

///
/// Solve $R x = \tilde{b}$.  Vector $b$ gets overwritten with $x$.
///
    const auto R = A.upper();
    blas::sv(NoTrans, R, _b);

    cout << "x = " << _b << endl;
}
