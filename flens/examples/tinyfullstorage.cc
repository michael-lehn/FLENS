#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int m = 3;
    const int n = 8;

    GeTinyMatrix<TinyFullStorage<double, m, n> >   A, B;


    A.engine().fill(1);
    A(2,3) = 10;

    A *= 2;
    cout << "A = " << A << endl;

    A /= 2;
    cout << "A = " << A << endl;

    B = A;
    cout << "B = " << B << endl;

    A = A + B;
    cout << "A = " << A << endl;
}
