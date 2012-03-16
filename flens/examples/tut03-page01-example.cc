#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;

    GEMatrix A(4,4);
    A = 1, 2, 3,  4,
        5, 6, 7,  8,
        9, 8, 7,  6,
        5, 4, 3, 20;

    cout << "A = " << A << endl;

    ///
    ///  ...
    ///
    auto U = A.upper();
    cout << "U = " << U << endl;

    ///
    ///  ...
    ///
    auto SU = U.symmetric();
    cout << "SU = " << SU << endl;

    ///
    ///  ...
    ///
    auto SL = A.lower().symmetric();
    cout << "SL = " << SL << endl;

    ///
    ///  ...
    ///
    auto L = A.lowerUnit();
    cout << "L = " << L << endl;

    return 0;
}
