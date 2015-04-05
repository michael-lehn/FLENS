#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4000;

    Matrix         A(n,n);
    Vector         x(n), b(n);
    IndexVector    piv(n);

    fillRandom(A);
    fillRandom(x);

    b = A*x;

    lapack::sv(A, piv, b);

    b -= x;

    cerr << "blas::asum(diff) = " << blas::asum(b) << endl;
}

