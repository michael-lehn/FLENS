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
/// Compute the factorization $A = LQ$.  Note that the workspace gets created
/// implicitly and temporarily.  So you might not want to do this inside a loop.
///
    lapack::lqf(A, tau);
    //lapack::lqf(A, tau, work);

///
/// Solve $L u = b$.  Vector $b$ gets overwritten with $u$.
///
    const auto L = A.lower();
    blas::sv(NoTrans, L, b);

///
/// Explicitly setup $Q$.
///
    Q = A;
    lapack::orglq(Q, tau);
    //lapack::orglq(Q, tau, work);

    cout << "Q = " << Q << endl;

///
/// Compute $x = Q^T u$.
///
    DenseVector<Array<double> >  x;
    x = transpose(Q)*b;

    cout << "x = " << x << endl;
}
