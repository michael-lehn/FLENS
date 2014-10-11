#include <cxxstd/iostream.h>
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

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;
    A *=I;

    b = 20,
       -33,
       -43,
        49;
    b *= I;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the factorization $A = QR$.
///
    lapack::qrf(A, tau);

///
/// Compute $\tilde{b} = Q^H b$. Vector $b$ gets overwritten with $\tilde{b}$.
///
    lapack::unmqr(Left, ConjTrans, A, tau, b);

///
/// Solve $R x = \tilde{b}$.  Vector $b$ gets overwritten with $x$.
///
    blas::sv(NoTrans, A.upper(), b);

    cout << "x = " << b << endl;
}
