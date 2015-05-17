#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >  DGeMatrix;

    const int n  = 9;

    DGeMatrix           A(n,n), B(n,n);

    fillRandom(A);

    cout << "A = " << A << endl;

    B = A*A;

    cout << "B = " << B << endl;

    B = A.lowerUnit()*B;

    cout << "B = " << B << endl;

    cout << "|| B ||_1 = " << blas::asum(B) << endl;
    return 0;
}
