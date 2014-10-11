#include <cxxstd/iostream.h>

///
/// Include our code snippet containing the high-level interface
///
#include "tut05-mylapack-version1.h"

using namespace flens;
using namespace std;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3.1;

    b = 20,
       -33,
       -43,
        49;

    cerr << "A = " << A << endl;
    cerr << "b = " << b << endl;

///
/// We now can use our interface just the way we are using __FLENS-LAPACK__ ...
///
/// :links: __FLENS-LAPACK__ -> doc:flens/lapack/lapack
///
    IndexType info = mylapack::trf(A, piv);

    if (info==0) {
///
///     ... for a user only the namspace differs.
///
        mylapack::trs(NoTrans, A, piv, b);

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }
}

