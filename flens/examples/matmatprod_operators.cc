#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4,4), B(4,4), C(4,4);

    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    B = 17, 18, 19, 20,
        21, 22, 23, 24,
        25, 26, 27, 28,
        29, 30, 31, 32;

    auto U = A.upper();
    auto S = A.upper().symmetric();

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "S = " << S << endl;
    cout << "U = " << U << endl;

//
//  compute  C = A*B
//
    C = A*B;

    cout << "C = A*B = " << C << endl;

//
//  compute  C = A^T*B
//
    C = transpose(A)*B;

    cout << "C = A^T*B = " << C << endl;

//
//  compute  C = S*B
//
    C = S*B;

    cout << "C = S*B = " << C << endl;

//
//  compute  B = U*B
//
    B = U*B;

    cout << "B = U*B = " << B << endl;
}
