//
// compile with:
//
// clang++ -std=c++0x qr_lapack.cc -I../../ -I/opt/local/include/ -L /opt/local/lib -lgmpxx -lgmp
//
//

#include <iostream>
#include <gmpxx.h>
#include <flens/flens.cxx>

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

using namespace std;
using namespace flens;

// typedef mpf_class   T;

typedef double   T;

#define TEST_CASE_1

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef GeMatrix::IndexType         IndexType;

    typedef DenseVector<Array<T> >      DenseVector;

    const Underscore<IndexType> _;

    const IndexType m = 4,
                    n = 5;

    GeMatrix            Ab(m, n);

#   ifdef TEST_CASE_1
    Ab =  3,  -1,   0,   0,   4,
         -1,   3,  -1,   0,  -2,
          0,  -1,   3,  -1,  -7,
          0,   0,  -1,   3,   8;
#   endif

#   ifdef TEST_CASE_2
    Ab =  2,   3,  -1,   0,  20,
         -6,  -5,   0,   2, -33,
          2,  -5,   6,  -6, -43,
          4,   6,   2,  -3,  49;
#   endif

#   ifdef TEST_CASE_3
    Ab =  2,  -6,   2,   4,  20,
          3,  -5,  -5,   6, -33,
         -1,   0,   6,   2, -43,
          0,   2,  -6,  -3,  49;
#   endif

    cerr << "Ab = " << Ab << endl;

    DenseVector tau(std::min(m,n));
    DenseVector work(n);

    lapack::qr2(Ab(_,_(1,m)), tau, work);
    cerr << "-> Ab = " << Ab << endl;

    lapack::qrs(Ab(_,_(1,m)), tau, Ab(_,_(m+1,n)), work);

    cerr << "-> Ab = " << Ab << endl;
}
