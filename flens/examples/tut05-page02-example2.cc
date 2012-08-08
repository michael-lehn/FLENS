#include <iostream>

#define USE_CXXLAPACK
#include <flens/flens.cxx>

///
/// We use our high-level interface that we defined before.
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

///
/// We define an underscore object to ease the creation of matrix/vector views.
///
    Underscore<IndexType>  _;

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

    IndexType info = mylapack::trf(A, piv);

    if (info==0) {
///
///     Passing `b(_(1,n))` as last argument to `trs` will cause a compile-time
///     error.  This is because the vector view `b(_(1,n))` is a temporary
///     object but `trs` expects a lvalue reference.
///
        mylapack::trs(NoTrans, A, piv, b(_(1,n)));
//
//      This would work (but it might be inconvenient):
//
//      auto _b = b(_(1,n));
//      mylapack::trs(NoTrans, A, piv, _b);
//

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }
}

