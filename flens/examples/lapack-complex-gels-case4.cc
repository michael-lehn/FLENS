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

    GeMatrix     A(3, 4);
    DenseVector  b(4);

///
/// Setup the matrix $A$.
///
    A = 1, 4, 7, 10,
        2, 5, 8, 11,
        3, 6, 9, 13;

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
    b *= 4.0*I;

///
/// The least square solution $x$ will be stored in the first three components
/// of vector $b$.  So vectors view come in handy:
///
    auto x = b(_(1,3));

///
/// Solve $\min\limits_{x} \|A^H x - b \|$
///
    lapack::ls(ConjTrans, A, b);

    cout << "x = " << x << endl;

    return 0;
}

