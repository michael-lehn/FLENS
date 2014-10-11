#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>                  Complex;
    const Complex                            I(0,1);

    typedef GeMatrix<FullStorage<Complex> >  GeMatrix;
    typedef DenseVector<Array<Complex> >     DenseVector;
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
/// Make the matrix/vectors somehow complex
///
    A *= I;
    b *= 2.0*I;

///
/// Find the minimal norm solution of $A^H x = b$.
///
    lapack::ls(ConjTrans, A, x);

    cout << "x = " << x << endl;

    return 0;
}

