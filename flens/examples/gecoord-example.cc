#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
///
/// Define a matrix with coordinate storage
///
    GeCoordMatrix<CoordStorage<double> >    A(5, 5);

///
/// When assembling a stiffness matrix we accumulate values.  Note that
/// an assignmen like 'A(2,3) = 4' is not allowed.
///
    A(3,4) += 4;
    A(3,4) += 4;
    A(2,4) += 4;
    A(3,4) += 4;

///
/// You see that values are justed added to a list.
///
    cout << "A = " << A << endl;

///
/// We explicitly can trigger the actual accumulation where values at identical
/// coordinates get summed up.
///
    A.engine().accumulate();
    cout << "A = " << A << endl;

///
/// For debugging its sometime usefull to convert the sparse matrix to a
/// matrix with full storage.
///
    GeMatrix<FullStorage<double> >  B = A;
    cout << "B = " << B << endl;

///
/// Once the stiffness matrix is setup you usually want to convert the sparse
/// matrix from coordinate storage to a compress format, e.g. the 'compressed
/// row storage' (CRS) format:
///
    GeCRSMatrix<CRS<double> >  C = A;
    cout << "C = " << C << endl;

///
/// Of course we also can convert the CRS matrix to GeMatrix for debugging.
///
    GeMatrix<FullStorage<double> >  D = C;
    cout << "D = " << D << endl;

///
/// Once the matrix is converted to the CRS format we can compute matrix-vector
/// products efficiently.
///
    DenseVector<Array<double> > x(5), y;

    x = 1;
    //y = D*x;
    blas::mv(NoTrans, double(1), C, x, double(0), y);
    //y = D*x;

    cout << "y = " << y << endl;
}
