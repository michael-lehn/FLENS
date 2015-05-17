#include <flens/flens.cxx>
#include <iostream>

using namespace std;
using namespace flens;

int
main()
{
///
/// Define a $3 \times 3$ triangular matrix with memory for $9$ doubles.
///
    TrMatrix<FullStorage<double> >    L(3, Lower);

///
/// Init L somehow
///
    L(1,1) = 11;
    L(2,1) = 21; L(2,2) = 22;
    L(3,1) = 31; L(3,2) = 32; L(3,3) = 33;

///
/// The following line is equivalent to:
/// `TrMatrix<FullStorage<double> >::GeneralView  A = L.general();`
///
    auto  A = L.general();

///
/// Change any element from the upper or lower triangular part
///
    A(1,3) = 42;
    A(2,2) = 42*2;
    A(3,1) = 42*3;

///
/// Note the output operator for `TrMatrix` only reads elements from the lower
/// triangular part from memory.  Elements from the upper part are just assumed
/// to be zero.
///
    cout << "A = " << A << endl;
    cout << "L = " << L << endl;
}
