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
    typedef GeMatrix<FullStorage<double, ColMajor> >      DGeMatrix;
    typedef GeMatrix<FullStorageView<double, ColMajor> >  DGeMatrixView;

///
/// We define an underscore operator for selecting index ranges.  The template
/// class `Underscore` requires the index type of `DGeMatrix` as template
/// parameter.
///
    const Underscore<DGeMatrix::IndexType>  _;

///
/// Then we setup some matrix.
///
    DGeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

///
/// Create a view `B` that references from `A` the second and third row.
///
    DGeMatrixView B = A(_(2,3),_);

    cout << "B = " << B << endl;

///
/// We change an element of `B`.  So we also change `A`.
///
    B(1,2) = 42;

    cout << "Changed B(1,2)." << std::endl;
    cout << "A = " << A << endl;


///
/// Create a copy `C` of `B`
///
    DGeMatrix C = B;

    cout << "C is a copy of B." << std::endl;
    cout << "C = " << C << endl;

///
/// We change an element of `C`.  So this neither will change `A` nor `B` as
/// `C` has its own memory.
///
    C(1,2) = 666;

    cout << "Changed C(1,2)." << std::endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
}
