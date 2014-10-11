#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >   GeMatrix;
    typedef DenseVector<Array<double> >      DenseVector;
    typedef typename DenseVector::IndexType  IndexType;

    const Underscore<IndexType>  _;

    GeMatrix     A(4, 3);
    DenseVector  x(4);

///
/// Setup the matrix $A$.
///
    A = 1, 5,  9,
        2, 6, 10,
        3, 7, 11,
        4, 8, 12;

///
/// Initially the first the elements of $x$ contain the right hand side vector
/// $b$.  For convenience we define a vector view $b$ that references these
/// elements.
///
    auto b = x(_(1,3));

///
/// Setup vector $b$.
///
    b = 30,
        70,
       110;

///
/// Find the minimal norm solution of $A^T x = b$.
///
    lapack::ls(Trans, A, x);

    cout << "x = " << x << endl;

    return 0;
}

