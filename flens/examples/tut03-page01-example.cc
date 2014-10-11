#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;

    GEMatrix A(4,4);
    A = 1, 2, 3,  4,
        5, 6, 7,  8,
        9, 8, 7,  6,
        5, 4, 3, 20;

    cout << "A = " << A << endl;

    ///
    /// `U` is a triangular matrix referencing the upper triangular part.
    ///
    auto U = A.upper();
    cout << "U = " << U << endl;

    ///
    /// `SU` is a symmetric matrix.  The upper triangular part of `SU` is
    /// referencing `U` (which in turn references the upper part of `A`).
    ///
    auto SU = U.symmetric();
    cout << "SU = " << SU << endl;

    ///
    /// `SL` is a symmetric matrix.  The lower triangular part of `SL` is
    /// referencing the lower part of `A`).
    ///
    auto SL = A.lower().symmetric();
    cout << "SL = " << SL << endl;

    ///
    /// `L` is a unit triangular matrix referencing the lower part of `A`.
    /// The term "unit" denotes that elements on the diagonal are assumed to
    /// be one.
    ///
    auto L = A.lowerUnit();
    cout << "L = " << L << endl;

    return 0;
}
