#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4,4);
    DenseVector<Array<double> >      x(4), y(4);

    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    x =  1,  2,  3,  4;

    cout << "A = " << A << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

//
//  compute  y = A*x
//
    blas::mv(NoTrans, 1.0, A, x, 0.0, y);

    cout << "y = A*x = " << y << endl;

//
//  compute  y = A^T*x
//
    blas::mv(Trans, 1.0, A, x, 0.0, y);

    cout << "y = A^T*x = " << y << endl;

//
//  compute  y = 1.5*y + A^T*x
//
    blas::mv(Trans, 1.0, A, x, 1.5, y);

    cout << "y = 1.5*y + A^T*x = " << y << endl;
}
