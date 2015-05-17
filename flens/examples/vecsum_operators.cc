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
    v = x + y + z;

    cout << "v = x + y + z = " << v << endl;

//
//  compute v = 2*x + 3*y + 4*z
//
    v = 2*x + 3*y + 4*z;

    cout << "v = 2*x + 3*y + 4*z = " << v << endl;
}
