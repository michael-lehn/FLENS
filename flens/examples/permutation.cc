#include <cmath>
#include <iostream>
#include <flens/examples/permutation.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef DenseVector<Array<double> >     RealDenseVector;
    typedef GeMatrix<FullStorage<double> >  RealGeMatrix;

    const int   n = 5;

    Permutation P(n);
    P.vector = 5, 3, 2, 1, 4;

    RealDenseVector                 x(5), y;

    x = 1, 2, 3, 4, 5;

    y = P*x;

    cout << "x = " << x << endl;
    cout << "y = P*x = " << y << endl;

    x = transpose(P)*y;
    cout << "x = P^T*y = " << x << endl;


    RealGeMatrix  I(n, n), B(n,n);
    I = 0;
    I.diag(0) = 1;

    B = P*I;
    cout << "B = P*I = " << B << endl;
}
