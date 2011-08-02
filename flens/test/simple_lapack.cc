// clang++ simple_lapack.cc -I../../ -I/opt/local/include/ -L /opt/local/lib -lgmpxx -lgmp

#include <iostream>
#include <gmpxx.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

// typedef mpf_class   T;

typedef mpq_class   T;

// typedef double   T;

#define TEST_CASE_2

int
main()
{
    typedef GeMatrix<FullStorage<T> >       GeMatrix;
    typedef GeMatrix::IndexType             IndexType;

    typedef DenseVector<Array<IndexType> >  IndexVector;

    const Underscore<IndexType> _;

    const IndexType m = 4,
                    n = 5;

    GeMatrix        Ab(m, n);
    IndexVector     iPiv(m);

#   ifdef TEST_CASE_1
    Ab =  3, -1,  0,  0,  4,
        -1,  3, -1,  0, -2,
         0, -1,  3, -1, -7,
         0,  0, -1,  3,  8;
#   endif

#   ifdef TEST_CASE_2
    Ab =  2,   3,  -1,   0,  20,
        -6,  -5,   0,   2, -33,
         2,  -5,   6,  -6, -43,
         4,   6,   2,  -3,  49;
#   endif

    cerr << "Ab = " << Ab << endl;

    /*
    auto A = Ab(_,_(1,m));
    auto b = Ab(_,_(m+1,n));

    lapack::sv(A, iPiv, b);
    */

    lapack::sv(Ab(_,_(1,m)), iPiv, Ab(_,_(m+1,n)));

    cerr << "-> Ab = " << Ab << endl;
    cerr << "-> iPiv = " << iPiv << endl;
}
