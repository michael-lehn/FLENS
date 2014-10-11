#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef double                               T;
    typedef DenseVector<Array<T> >               DEVector;
    typedef GeMatrix<FullStorage<T, ColMajor> >  GEMatrix;

    const T  alpha = 1.5,
             beta = 2.5;

    DEVector x(3), y(3), z(3);
    x = 1, 2, 3;
    y = 2, 3, 4;
    z = 3, 4, 5;

    GEMatrix A(3,3);
    A = 1, 2, 3,
        5, 6, 7,
        5, 4, 3;

    ///
    /// Compute $y = \beta y + \alpha A^T x$
    ///
    blas::mv(Trans, alpha, A, x, beta, y);

    ///
    /// Compute the update $y = y + z$
    ///
    blas::axpy(T(1), z, y);

    cout << "y = " << y << endl;

    return 0;
}
