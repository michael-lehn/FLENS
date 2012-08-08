#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    typedef complex<double>             Complex;
    const Complex                       I(0,1);


    GeMatrix<FullStorage<Complex> >     A(4,4);
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
/// Compute $\tilde{b} = Q^T b$.  Note that we do not setup $A$ explicitly.
///
    lapack::unmqr(Left, ConjTrans, A, tau, b);
    // lapack::ungqr(A, tau, b, work);

///
/// Solve $R x = \tilde{b}$.  Vector $b$ gets overwritten with $x$.
///
    const auto R = A.upper();
    blas::sv(NoTrans, R, b);

    cout << "x = " << b << endl;
}
