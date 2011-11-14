//
// compile with:
//
// clang++ -std=c++0x simple_vec.cc -I../../
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
    typedef DenseVector::IndexType          IndexType;

    const Underscore<IndexType>             _;

    DenseVector x(5), y(5);
    x = 1, 2, 3, 4, 5;

    cout << "x =      " << x << endl;

    y = x.reverse();
    cout << "y =      " << y << endl;


    cout << "x(_(2,5)) = " << x(_(2,5)) << endl;

    auto view       = x(_(2,5));
    auto viewView   = view(_(3,4));

    view = 9;
    viewView = 4;
    cout << "x =      " << x << endl;
}