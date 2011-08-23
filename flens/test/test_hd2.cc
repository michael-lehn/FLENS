//
// compile with:
//
// clang++ -std=c++0x test_hd2.cc -I../../
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

#define TEST_CASE_1

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef GeMatrix::IndexType         IndexType;

    typedef DenseVector<Array<T> >      DenseVector;

    const Underscore<IndexType> _;

    const IndexType n = 4;

    GeMatrix            A(n, n), Z(n, n);
    DenseVector         wr(n), wi(n);

#   ifdef TEST_CASE_1
    A =  3,  -1,   0,   0,
        -1,   3,  -1,   0,
         0,  -1,   3,  -1,
         0,   0,  -1,   3;
#   endif

#   ifdef TEST_CASE_2
    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;
#   endif

#   ifdef TEST_CASE_3
    A =  2,  -6,   2,   4,
         3,  -5,  -5,   6,
        -1,   0,   6,   2,
         0,   2,  -6,  -3;
#   endif

    cerr << "A = " << A << endl;

    DenseVector tau(n-1);
    DenseVector work(n);

    lapack::hd2(IndexType(1), n, A, tau, work);

    cerr << "A = " << A << endl;

    IndexType info = lapack::hseqr(lapack::HSEQR::Eigenvalues, lapack::HSEQR::Init,
                                   IndexType(1), n, A, wr, wi, Z, work);

    cerr << "-> info =  " << info << endl;
    cerr << "-> A =  " << A << endl;
    cerr << "-> wr = " << wr << endl;
    cerr << "-> wi = " << wi << endl;
    cerr << "-> Z =  " << Z << endl;
}
