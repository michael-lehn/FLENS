#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;
    GEMatrix A(4,4);
    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    cout << "A = " << A << endl;

    ///
    ///  We define a underscore object such that `_(from,to)` later defines
    ///  a index range.
    ///
    Underscore<GEMatrix::IndexType> _;

    ///
    ///  Create a view of a submatrix and print it.
    ///
    cout << "A(2:3,1:3) = " << A(_(2,3),_(1,3)) << endl;

    ///
    ///  Create a view and give it a name.  Note that `A11` has type
    ///  `GeMatrix<FullStorageView<double, ColMajor> >`.  You see, the `auto`
    ///  type makes life much easier here.
    ///
    auto A11 = A(_(1,2),_(1,2));

    ///
    ///  Instead of `auto` we also could have used the internal `View` typedef
    ///  of `GeMatrix`:
    ///
    GEMatrix::View A22 = A(_(3,4),_(3,4));

    ///
    ///  Copy A22 to A11
    ///
    A11 = A22;

    ///
    /// Copy unnamed sub-blocks of A (from top-right to bottom-left).
    ///
    A(_(3,4),_(1,2)) = A(_(1,2),_(3,4));

    ///
    ///  As `A11` and `A22` are just referencing parts of `A` we actually
    ///  changed `A`
    ///

    cout << "A = " << A << endl;
    return 0;
}
