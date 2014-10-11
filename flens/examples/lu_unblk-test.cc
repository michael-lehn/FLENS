#include <cxxstd/algorithm.h>
#include <flens/flens.cxx>
#include <flens/examples/lu_unblk.h>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >   RealGeMatrix;

    const int m = 2000;
    const int n = 2000;
    const int mn = std::min(m, n);

///
/// Setup test case and make a copy of the input
///
    RealGeMatrix  A(m, n), A_check;
    fillRandom(A);
    A_check = A;

///
///  Factorize $A$.
///
    if (lu_unblk(A)) {
        cerr << "A is singular" << endl;
        return 1;
    }

///
/// Compute $|A - L*U|_1$ by brute force.  Note that $A$ is not required to
/// be square.
///
    Underscore<int> _;
    RealGeMatrix L = A(_,_(1,mn)).lowerUnit();
    RealGeMatrix U = A(_(1,mn),_).upper();

    A_check -= L*U;

    cout << "asum(A-L*U) = " << blas::asum(A_check) << endl;
    return 0;
}
