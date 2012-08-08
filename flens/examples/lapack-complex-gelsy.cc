#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>                  Complex;
    const Complex                            I(0,1);

    typedef GeMatrix<FullStorage<Complex> >  GeMatrix;
    typedef typename GeMatrix::IndexType     IndexType;
    typedef DenseVector<Array<IndexType> >   IndexVector;
    typedef DenseVector<Array<Complex> >     DenseVector;

    const Underscore<IndexType>  _;

    GeMatrix     A(3, 4);
    DenseVector  x(4);
    IndexVector  jPiv(4);

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
/// Make $A$ and $b$ somehow complex ;-)
///
    A *= I;
    b *= 3.0*I;

///
/// All columns are free columns.
///
    jPiv = 0;

///
/// We want that the condition of the leading submatrix $R_{11}$ is less than
/// `1/rCond`.
///
    double    rCond = 1E-8;

///
/// Find the minimal norm solution of $Ax = b$.
///
    IndexType rank = lapack::lsy(A, x, jPiv, rCond);

    cout << "x = " << x << endl;

///
/// Output of permuation $P$ and the effective rank.
///
    cout << endl << endl << endl
         << "Some additional information:"
         << endl;

    cout << "jPiv = " << jPiv << endl;
    cout << "rank = " << rank << endl;

    return 0;
}

