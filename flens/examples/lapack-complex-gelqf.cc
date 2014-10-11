#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

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
    A *=I;

    b = 20,
       -33,
       -43,
        49;
    b *= I;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the factorization $A = LQ$.  Note that the workspace gets created
/// implicitly and temporarily.  So you might not want to do this inside a loop.
///
    lapack::lqf(A, tau);
    // lapack::lqf(A, tau, work);

///
/// Solve $L u = b$.  Vector $b$ gets overwritten with $u$.
///
    blas::sv(NoTrans, A.lower(), b);

///
/// Compute $x = Q^T u$. Vector $b$ gets overwritten with $x$.
///
    lapack::unmlq(Left, ConjTrans, A, tau, b);
    // lapack::unmlq(Left, ConjTrans, A, tau, b, work);

    cout << "x = " << b << endl;
}
