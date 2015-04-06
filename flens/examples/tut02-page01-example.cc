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
    cxxblas::gemv(A.numRows(), A.numCols(),
                  alpha,
                  true, false, A.data(), A.strideRow(), A.strideCol(),
                  x.data(), x.stride(),
                  beta,
                  y.data(), y.stride());

    ///
    /// Compute the update $y = y + z$
    ///
    cxxblas::axpy(y.length(), T(1), z.data(), z.stride(), y.data(), y.stride());

    cout << "y = " << y << endl;

    return 0;
}
