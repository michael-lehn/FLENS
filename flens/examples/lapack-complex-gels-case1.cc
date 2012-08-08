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
    typedef DenseVector<Array<Complex> >     DenseVector;
    typedef typename DenseVector::IndexType  IndexType;

    const Underscore<IndexType>  _;

    GeMatrix     A(4, 3);
    DenseVector  b(4);

///
/// Setup the matrix $A$.
///
    A =  1,  2,  3,
         4,  5,  6,
         7,  8,  9,
        10, 11, 13;

///
/// Setup vector $b$
///
    b =  6,
        15,
        24,
        34;

///
/// Make the matrix/vectors somehow complex
///
    A *= I;
    b *= I;

///
/// The least square solution $x$ will be stored in the first three components
/// of vector $b$.  So vectors view come in handy:
///
    auto x = b(_(1,3));

///
/// Solve $\min\limits_{x} \|Ax - b \|$
///
    lapack::ls(NoTrans, A, b);

    cout << "x = " << x << endl;

    return 0;
}

