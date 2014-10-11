#include <cxxstd/iostream.h>

#define USE_PLAYGROUND
#include <flens/flens.cxx>

using namespace std;
using namespace flens;
using namespace dft;

typedef complex<double>   T;

int
main(int argc, char* argv[])
{
    ///
    ///  Define some convenient typedefs for the vector types
    ///
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;
    Matrix A(n, n);
    IndexVector Pivots;

    ///
    /// Fill in random values
    ///
    fillRandom(A);
    cout << " A = " << A << endl;

    ///
    /// Calculate determinant based on LU factorization
    ///
    cout << lapack::extensions::det(A, Pivots) << endl;

    return 0;
}
