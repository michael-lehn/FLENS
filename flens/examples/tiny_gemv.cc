#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef double D;
    typedef GeTinyMatrix<TinyFullStorage<D, 4, 4> >  TMatrix_4_4;
    typedef TinyVector<TinyArray< double, 4 > >      TVector_4;

    TMatrix_4_4 A;
    A(1,1) = 1;  A(1,2) = 2;  A(1,3) = 3;  A(1,4) =  4;
    A(2,1) = 5;  A(2,2) = 6;  A(2,3) = 7;  A(2,4) =  8;
    A(3,1) = 9;  A(3,2) = 8;  A(3,3) = 7;  A(3,4) =  6;
    A(4,1) = 5;  A(4,2) = 4;  A(4,3) = 3;  A(4,4) = 20;

    cout << "A = " << A << endl;

    TVector_4 v;
    v.fill(1);

    cout << "v = " << v << endl;

    cout << "TVector_4: v = " << v << endl;
    cout << "sizeof(v)= " << sizeof(v) << endl;

    TVector_4  w;

    w = A * v;
    cout << "w = A*v = " << w << endl;

    w = transpose(A) * v;
    cout << "w = A*v = " << w << endl;
}
