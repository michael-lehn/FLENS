//
// compile with:
//
// clang++ -std=c++0x simple_matvec.cc -I../../ -I/opt/local/include/ -L /opt/local/lib -lgmpxx -lgmp
//
//

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

    GeMatrix A(3,3);
    A = 1, 2, 3,
        3, 2, 1,
        1, 2, 1;

    cout << "A = " << A << endl;
    auto _x = A(_(2,3),2);
    cout << "_x = A(2:3,2) = " << _x << endl;
    cout << "_x.firstIndex() = " << _x .firstIndex() << endl;
}