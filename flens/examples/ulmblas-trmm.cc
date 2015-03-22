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

    GeMatrix<FullStorage<complex<float>, RowMajor> >   A(m,m);
    GeMatrix<FullStorage<complex<double>, ColMajor> >  B(m,n);
    GeMatrix<FullStorage<complex<double>, ColMajor> >  C(n,m);

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

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

    B = A.upper()*B;
    cout << "B = " << B << endl;

    B = A.lowerUnit()*B;
    cout << "B = " << B << endl;

    C = C * A.upper();
    cout << "C = " << C << endl;
}
