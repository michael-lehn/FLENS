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

    GeMatrix     A(3, 4);
    DenseVector  x(4);

///
/// Setup the $3 \times 4$ matrix $A$.
///
    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;

///
/// The right hand side vector $b$ will be stored in vector $x$.
///
    auto b = x(_(1,3));
    b =  30,
         70,
        110;

///
/// Find the minimal norm solution of $Ax = b$.
///
    lapack::ls(NoTrans, A, x);

    cout << "x = " << x << endl;

    return 0;
}

