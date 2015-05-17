#include <flens/flens.cxx>
#include <algorithm>
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
    cxxblas::copy(v.length(), x.data(), x.stride(), v.data(), v.stride());

    // v += y
    cxxblas::axpy(v.length(), 1., y.data(), y.stride(), v.data(), v.stride());

    // v += z
    cxxblas::axpy(v.length(), 1., z.data(), z.stride(), v.data(), v.stride());

    cout << "v = x + y + z = " << v << endl;

//
//  compute v = 2*x + 3*y + 4*z
//

    // set all element of v to zero
    std::fill_n(v.data(), v.length(), 0.0);

    // v += 2*x
    cxxblas::axpy(v.length(), 2., x.data(), x.stride(), v.data(), v.stride());

    // v += 3*y
    cxxblas::axpy(v.length(), 3., y.data(), y.stride(), v.data(), v.stride());

    // v += 4*z
    cxxblas::axpy(v.length(), 4., z.data(), z.stride(), v.data(), v.stride());

    cout << "v = 2*x + 3*y + 4*z = " << v << endl;
}
