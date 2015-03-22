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
    const int k = 5;
    const int n = 3;

    const complex<float> x = double(1);

    //typedef typename std::common_type<double, complex<float> >::type   T;
    //typedef typename std::common_type<complex<float>, double>::type   T;
    //typedef typename std::common_type<float, double>::type   T;

    //foo(T(0));

    GeMatrix<FullStorage<complex<float>, RowMajor> >  A(m,k);
    GeMatrix<FullStorage<double, ColMajor> >  B(k,n);
    GeMatrix<FullStorage<complex<double>, RowMajor> >  C(m,n);
    GeMatrix<FullStorage<complex<double>, RowMajor> >  C2;

    int count = 0;

    // init A
    for (int i=1; i<=m; ++i) {
        for (int j=1; j<=k; ++j) {
            A(i,j) = ++count;
        }
    }
    // init B
    for (int i=1; i<=k; ++i) {
        for (int j=1; j<=n; ++j) {
            B(i,j) = ++count;
        }
    }
    // init C
    for (int i=1; i<=m; ++i) {
        for (int j=1; j<=n; ++j) {
            C(i,j) = ++count;
        }
    }

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

    C = A*B;

    cout << "C = " << C << endl;

    C2 = transpose(B)*conjTrans(A);

    cout << "C2 = " << C2 << endl;
}
