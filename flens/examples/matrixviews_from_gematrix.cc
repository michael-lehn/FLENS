#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(3,4);

    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;

    cout << "A = " << A << endl;
    cout << "A.upper() = " << A.upper() << endl;
    cout << "A.upperUnit() = " << A.upperUnit() << endl;
    cout << "A.lower() = " << A.lower() << endl;
    cout << "A.lowerUnit() = " << A.lowerUnit() << endl;
}
