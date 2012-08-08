#include <iostream>

#define USE_CXXLAPACK
#include <flens/flens.cxx>

// We now include the final version of `mylapack`.
#include "tut05-mylapack-version3.h"

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

    const Underscore<IndexType> _;

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

    IndexType info = mylapack::trf(A, piv);

    if (info==0) {
        const auto b_old = b;

//      We call trs with a lvalue ...
        mylapack::trs(NoTrans, A, piv, b);

        b = b_old;

//      ... and a rvalue. Just to see if now both cases work.
        mylapack::trs(NoTrans, A, piv, b(_(1,n)));

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }
}

