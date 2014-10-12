#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

///
/// In the `main` function we again start with a typedef for a `GeMatrix` with
/// elements of type `double` that are stored column-wise in memory.  Also,
/// we define a suitable range operator.
///
int
main()
{
    typedef FullStorage<double, ColMajor>  FS;
    typedef GeMatrix<FS>                   GeMatrix;
    const Underscore<GeMatrix::IndexType>  _;

///
/// Then we setup some matrix.
///
    GeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

///
/// Matrix `B` is of type `const GeMatrix::View`.  So in particular this matrix
/// is in const context.
///
    const auto B = A(_(1,3),_);

///
/// Matrix `C` is of type `GeMatrix::ConstView` but it is *not* in const
/// context.
///
    auto C = B(_(1,3),_);

///
/// Matrix `D` would be of type `GeMatrix::View` and would reference a part
/// of a matrix of type `GeMatrix::ConstView`.  This would violate *const
/// correctness*.
    auto D = C(_(1,3),_);

    // D references a const matrix view.  So this must trigger a compile
    // time error!
    D(1,1) = 666;
    cout << "D = " << D << endl;
}
