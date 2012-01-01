#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    ///
    ///  You can use the underscore range operator `_(from, to)` to set the
    ///  index range of a matrix.
    ///
    ///  Note that in previous versions of FLENS the underscore object `_` was
    ///  already defined as a global variable.  You now have to define it
    ///  manually.  This is because every matrix can now use  different
    ///  index types (usually int or long).
    ///
    typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;
    Underscore<GEMatrix::IndexType> _;
    GEMatrix A(_(0,3),_(-2,1));

    for (int i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (int j=A.firstCol(); j<=A.lastCol(); ++j) {
            A(i,j) = 100*i+j;
        }
    }

    ///
    ///  Or define a `GeMatrix` type with a different index base.  Here we set
    ///  the default base to 0 such that both, row and col indices, start at
    ///  zero.
    ///
    typedef IndexOptions<int, 0>                                 ZeroBased;
    typedef GeMatrix<FullStorage<double, ColMajor, ZeroBased> >  GEMatrixZB;
    GEMatrixZB  B(3,3);

    for (int i=B.firstRow(); i<=B.lastRow(); ++i) {
        for (int j=B.firstCol(); j<=B.lastCol(); ++j) {
            B(i,j) = 100*i+j;
        }
    }


    ///
    ///  Print matrix dimensions and content.  Note that this time we use
    ///  methods `rows()` and `cols()` which return range objects.
    ///
    cout << "A(" << A.rows() << ", " << A.cols() << ") = " << A << endl;
    cout << "B(" << B.rows() << ", " << B.cols() << ") = " << B << endl;

    return 0;
}
