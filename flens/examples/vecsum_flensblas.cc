#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    DenseVector<Array<double> >  x(4), y(4), z(4), v(4);

    x =  1,  2,  3,  4;
    y =  5,  6,  7,  8;
    z =  9, 10, 11, 12;

    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

//
//  compute  v = x + y + z
//

    // v = x
    blas::copy(x, v);

    // v += y
    blas::axpy(1., y, v);

    // v += z
    blas::axpy(1., z, v);

    cout << "v = x + y + z = " << v << endl;

//
//  compute v = 2*x + 3*y + 4*z
//

    // set all element of v to zero
    v = 0.0;

    // v += 2*x
    blas::axpy(2., x, v);

    // v += 3*y
    blas::axpy(3., y, v);

    // v += 4*z
    blas::axpy(4., z, v);

    cout << "v = 2*x + 3*y + 4*z = " << v << endl;
}
