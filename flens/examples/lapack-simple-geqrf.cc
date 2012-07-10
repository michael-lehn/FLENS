#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    GeMatrix<FullStorage<double> >     A(4,4);
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
/// Compute $\tilde{b} = Q^H b$. Vector $b$ gets overwritten with $\tilde{b}$.
///
    lapack::ormqr(Left, Trans, A, tau, b);
    //lapack::ormqr(Left, Trans, A, tau, b, work);

///
/// Solve $R x = \tilde{b}$.  Vector $b$ gets overwritten with $x$.
///
    blas::sv(NoTrans, A.upper(), b);

    cout << "x = " << b << endl;
}
