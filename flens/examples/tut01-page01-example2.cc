#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

///
/// In the `main` function we again start with a typedef for a `GeMatrix` with
/// elements of type `double` that are stored column-wise in memory.  Two other
/// typedefs are used for corresponding matrix views and matrix const-views.
/// For nicely selecting rectangular parts we define a suitable range operator.
///
int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >      DGeMatrix;
    typedef DGeMatrix::View                               DGeMatrixView;
    typedef DGeMatrix::ConstView                          DGeMatrixConstView;

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
/// Create a view `B` that references from `A` the second and third row. And
/// a const-view `C` referencing the second and third column of `B`.
///
    DGeMatrixView             B = A(_(2,3),_);
    const DGeMatrixConstView  C = B(_,_(2,3));

    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

///
/// Changing elements through a const-view gets prevented by the compiler. So
/// the following line of code would cause a compile-time error:
///
    // C(1,1) = 31415;  // error

///
/// You can create a view or const-view from a view   But from a const-view
/// you only can create a const-view.  Attempting to create a non-const view
/// from a const view will give a compile-time error:
///
    // DGeMatrixView   X = C(_(1,2),_(1,2));  // error

///
/// A `GeMatrix` instance in const-context returns a const-view.  As you can not
/// construct a view from a const-view this would cause a compile-time error
///
    // const DGeMatrix  Z(3,4);
    // DGeMatrixView    B = A(_(2,3),_);  // erro

///
/// We change an element of `B`.  So we also change `A` and `C`.
///
    B(1,2) = 42;

    cout << "Changed B(1,2)." << std::endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;


///
/// Create a copy `D` of `B`
///
    DGeMatrix D = B;

    cout << "D is a copy of B." << std::endl;
    cout << "D = " << D << endl;

///
/// We change an element of `D`.  So this neither will change `A` nor `B` as
/// `D` has its own memory.
///
    D(1,2) = 666;

    cout << "Changed D(1,2)." << std::endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
    cout << "D = " << D << endl;
}
