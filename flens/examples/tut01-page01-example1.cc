///
///  We simply include everything of FLENS
///
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    ///
    ///  Typedef for a general matrix with elements of type `double`.
    ///  Internally elements are stored in *column major* order.
    ///
    typedef GeMatrix<FullStorage<double, ColMajor> >  DGeMatrix;

    ///
    ///  Matrix A gets dynamically allocated and then initialized.
    ///
    DGeMatrix A(4,4);
    A = 1, 2, 3,  4,
        5, 6, 7,  8,
        9, 8, 7,  6,
        5, 4, 3, 20;

    ///
    ///  Print the matrix content using output streams:
    ///
    cout << "A = " << A << endl;

    ///
    ///  We print some information about matrix dimensions and index ranges.
    ///  You will see that by default in FLENS indices start at 1 (like in
    ///  Fortran):
    ///
    cout << "Dim. of A: " << A.numRows() << " x " << A.numCols() << endl;
    cout << endl;
    cout << "Row indices: " << A.firstRow() << ".." << A.lastRow() << endl;
    cout << endl;
    cout << "Col indices: " << A.firstCol() << ".." << A.lastCol() << endl;
    cout << endl;

    ///
    /// Also for element access (write) we provide a Fortran-Style interface:
    ///
    A(3,2) = 42;

    ///
    /// The same for read access:
    ///
    cout << "changed element: A(3,2) = " << A(3,2) << endl;

    cout << endl;

    cout << "A = " << A << endl;

    ///
    /// You also can fill the whole matrix with a new value:
    ///
    A = 42;

    cout << "A = " << A << endl;

    return 0;
}
