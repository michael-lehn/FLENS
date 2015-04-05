#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef typename Matrix::IndexType          IndexType;

    Underscore<IndexType> _;

    int m = 10;
    int n = 12;

    Matrix         A(m,n);

    int count = 0;

    for (IndexType j=1; j<=n; ++j) {
        for (IndexType i=1; i<=m; ++i) {
            A(i,j) = ++count;
        }
    }

    cout << "A = " << A << endl;
    cout << A(_(2,2,5),_(3,2,7)) << endl;

    auto A_ = A(_(1,2,5),_(3,2,7));

    Matrix B(A_.numRows(), A_.numCols());

    B = A_;

    cout << "B = " << B << endl;
}

