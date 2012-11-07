#include <algorithm>
#include <flens/flens.cxx>
#include <flens/examples/lu_blk.h>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >   RealGeMatrix;

    const int m = 2000;
    const int n = 2000;
    const int mn = std::min(m, n);

    RealGeMatrix  A(m, n), A_check;

    fillRandom(A);
    A_check = A;


    if (lu_blk(A)) {
        cerr << "A is singular" << endl;
        return 1;
    }

    Underscore<int> _;
    RealGeMatrix L = A(_,_(1,mn)).lowerUnit();
    RealGeMatrix U = A(_(1,mn),_).upper();

    A_check -= L*U;

    cout << "asum(A-L*U) = " << blas::asum(A_check) << endl;
    return 0;
}
