#include <flens/flens.cxx>
#include <flens/examples/lu/tile.h>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >   DGeMatrix;
    typedef DenseVector<Array<double> >      DDenseVector;

    const int m  = 10;
    const int n  = 10;
    const int bs = 6;

    DGeMatrix       A(m, n);

    int count = 1;
    for (int j=1; j<=n; ++j) {
        for (int i=1; i<=m; ++i) {
            A(i,j) = count++;
        }
    }

    cout << "A = " << A << endl;

    DDenseVector            buffer;
    TiledCopy<DGeMatrix>    A_(A, bs, buffer);

    cout << "Tiled Matrix A_" << endl;
    for (int i=1; i<=A_.numTileRows(); ++i) {
        for (int j=1; j<=A_.numTileCols(); ++j) {
            cout << "(" << i << ", " << j << "):" << A_(i,j) << endl;
        }
    }

    DGeMatrix B(m,n);

    A_.untile(B);

    cout << "B = " << B << endl;
}
