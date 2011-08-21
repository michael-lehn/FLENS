//
// compile with:
//
// clang++ -std=c++0x test_lanv2.cc -I../../
//
//

#include <iostream>
#include <flens/flens.cxx>

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

using namespace std;
using namespace flens;

// typedef mpf_class   T;

typedef double   T;

#define TEST_CASE_2

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef GeMatrix::IndexType         IndexType;

    typedef DenseVector<Array<T> >      DenseVector;

    const Underscore<IndexType> _;

    const IndexType n = 2;

    GeMatrix            A(n, n), V(n, n);
    T                   cs, sn;

#   ifdef TEST_CASE_1
    A =  3,  -1,
        -1,   3;
#   endif

#   ifdef TEST_CASE_2
    A =  2,   3,
        -6,  -5;
#   endif

#   ifdef TEST_CASE_3
    A =  2,  -1,
         1,   2;
#   endif

#   ifdef TEST_CASE_4
    A =  2,  -1,
        -1,   2;
#   endif

    cerr << "A = " << A << endl;

    lapack::lanv2(A(1,1), A(1,2), A(2,1), A(2,2),
                  V(1,1), V(1,2), V(2,1), V(2,2),
                  cs, sn);
    
    cerr << "-> A = " << A << endl;
    cerr << "-> V = " << V << endl;
    cerr << "-> cs = " << cs << ", sn = " << sn << std::endl;
}
