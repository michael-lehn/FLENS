#if defined(__SSE3__)
//#   undef __SSE3__
#endif


#include <cxxstd/iostream.h>
///
///  With header __flens.cxx__ all of FLENS gets included.
///
///  :links:  __flens.cxx__ -> file:flens/flens.cxx



#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    ///
    ///  Define some convenient typedefs for the matrix/vector types
    ///  of our system of linear equations.
    ///
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    ///
    ///  We also need an extra vector type for the pivots.  The type of the
    ///  pivots is taken for the system matrix.
    ///
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    ///
    ///  Set up the  problem ...
    ///
    const IndexType n = 4000;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    fillRandom(A);
    fillRandom(b);


    ///
    /// And solve it with __lapack::sv__
    ///
    /// :links: __lapack::sv__ -> doc:flens/lapack/ge/sv
    lapack::sv(A, piv, b);

}

