#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
///
/// We start with a typedef for a `GeMatrix` with elements of type `double`
/// that are stored column-wise in memory.
///
    typedef FullStorage<double, ColMajor>  FS;
    typedef GeMatrix<FS>                   GeMatrix;

///
/// Then we define matrix view types for referencing matrix parts.
///
    typedef GeMatrix::ConstView            GeMatrixConstView;
    typedef GeMatrix::View                 GeMatrixView;

///
/// We setup some matrix.
///
    GeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

///
/// We also define a "range-operator" to ease the selection of matrix parts.
/// for this purpose we simply define an object called "`_`" of type
/// `Underscore<I>` where `I` denotes the type of index used in the matrix.
///
    const Underscore<GeMatrix::IndexType>   _;

///
/// Now we let `B` reference the matrix part `A(2:3,2:4)`.  The range operator
/// lets you specify row/column ranges by `_(from, to)` similar to the Matlab
/// notation `from:to`.
///
    GeMatrixView       B = A(_(2,3),_(2,4));

///
/// Note that you can not define something like an uninitialized matrix view.
/// Just like you can not define uninitialized references in C++.  I.e. this
/// would give you an error
//
//  error:
//  GeMatrixView       B;
//

///
/// Next we define a matrix view that provides read-only access to `A(:,2,4)`.
/// So by the way, notice that you can also select all rows if you do not pass
/// arguments to the range operator "`_`".  Analogously you can select all
/// columns of course.
///
    GeMatrixConstView  C = A(_,_(2,4));

///
/// We overwrite `B` and by that `A`
///
    B = 111, 222, 333,
        444, 555, 666;

    cout << "Changed B:" << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

///
/// Note that `B` acts like a "regular" matrix and indices also start at One.
///
    B(1,1) = 999;

    cout << "Changed B(1,1):" << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

///
/// Of course you also can use unnamed matrix views:
///
    A(_,_(1,2)) = A(_,_(3,4));

    cout << "After copying A(:,1:2) = A(:,3:4)" << endl;
    cout << "A = " << A << endl;

///
/// And obviously it should not be possible to write to `C`.  Indeed this would
/// give you a compile-time error:
//
//  error:
//  C(1,1) = 3;
//  C = 1, 2,
//      3, 4,
//      5, 6;
}
