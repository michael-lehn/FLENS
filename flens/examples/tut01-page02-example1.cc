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
    ///  Typedef for a dense vector with elements of type `double`.
    ///
    typedef DenseVector<Array<double> >  DDenseVector;

    ///
    ///  Vector x gets dynamically allocated and then initialized.
    ///
    DDenseVector x(5);
    x = 1, 2, 3, 4, 5;

    ///
    ///  Print the vector content using output streams:
    ///
    cout << "x = " << x << endl;

    ///
    ///  We print some information about vector length and index ranges.
    ///  You will see that by default in FLENS indices start at 1 (like in
    ///  Fortran):
    ///
    cout << "Length of x: " << x.length() << endl;
    cout << endl;
    cout << "Indices: " << x.firstIndex() << ".." << x.lastIndex() << endl;
    cout << endl;

    ///
    /// Also for element access (write) we provide a Fortran-Style interface:
    ///
    x(2) = 42;

    ///
    /// The same for read access:
    ///
    cout << "changed element: x(2) = " << x(2) << endl;

    cout << endl;

    cout << "x = " << x << endl;

    ///
    /// You also can fill the whole vector with a new value:
    ///
    x = 42;

    cout << "x = " << x << endl;

    return 0;
}
