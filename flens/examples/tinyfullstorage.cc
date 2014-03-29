#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int m = 3;
    const int n = 8;

    GeTinyMatrix<TinyFullStorage<double, m, n> >   A, B;
    GeTinyMatrix<TinyFullStorageView<double, m, n/2, n> >   vA(A.data()+1);
    GeTinyMatrix<TinyConstFullStorageView<double, m, n/2, n> >   cvA(A.data());

    A.engine().fill(1);
    A(2,3) = 10;
        vA(1,1) = 4;
    cout << A << endl;
    cout << vA << endl;

    cout << cvA << endl;
    return 0;

    A *= 2;
    cout << "A = " << A << endl;

    A /= 2;
    cout << "A = " << A << endl;

    B = A;
    cout << "B = " << B << endl;

    A = A + B;
    cout << "A = " << A << endl;
}
