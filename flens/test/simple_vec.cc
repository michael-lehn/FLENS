#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef DenseVector<Array<double> >     DenseVector;

    DenseVector x(5), y(5);
    x = 1, 2, 3, 4, 5;

    cout << "x =      " << x << endl;

    y = x.reverse();
    cout << "y =      " << y << endl;
}