#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef DenseVector<Array<double> >     DenseVector;
    typedef GeMatrix<FullStorage<double> >  GeMatrix;

    Underscore<>  _;

    DenseVector x(5);
    x = 1, 2, 3, 4, 5;

    cout << "x =      " << x << endl;
    cout << "x(2:4) = " << x(_[2][4]) << endl;

    GeMatrix A(3,3);
    A = 1, 2, 3,
        3, 2, 1,
        1, 2, 1;

    cout << "A = " << A << endl;
}