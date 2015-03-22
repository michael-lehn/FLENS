#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

void
foo(int *i);

int
main()
{
    const int m = 4;
    const int n = 3;

    GeMatrix<FullStorage<float, RowMajor> >   A(m,m);
    GeMatrix<FullStorage<double, ColMajor> >  B(m,n);
    GeMatrix<FullStorage<double, ColMajor> >  C(m,n);
    GeMatrix<FullStorage<double, ColMajor> >  D(n,m), E;

    int count = 0;

    // init A
    for (int i=1; i<=A.numRows(); ++i) {
        for (int j=1; j<=A.numCols(); ++j) {
            A(i,j) = ++count;
        }
    }
    // init B
    for (int i=1; i<=B.numRows(); ++i) {
        for (int j=1; j<=B.numCols(); ++j) {
            B(i,j) = ++count;
        }
    }
    // init C
    for (int i=1; i<=C.numRows(); ++i) {
        for (int j=1; j<=C.numCols(); ++j) {
            C(i,j) = ++count;
        }
    }

    D = transpose(C);

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

    C = A.upper().symmetric()*B;
    cout << "C = " << C << endl;

    cout << "D = " << D << endl;
    cout << "A.upper().symmetric() = " << A.upper().symmetric() << endl;

    E = D * A.upper().symmetric();
    cout << "E = " << E << endl;
}
