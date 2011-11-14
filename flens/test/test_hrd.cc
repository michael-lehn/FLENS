//
// compile with:
//
// clang++ -std=c++0x -DWITH_REFBLAS -c -DCHECK_CXXLAPACK test_hrd.cc  -I../../
// gfortran test_hrd.o -lstdc++ tmp/*.o drivers.o ../lapack/interface/cblas_REF.a ../lapack/interface/blas_REF.a

#include <iostream>
#define FLENS_DEFAULT_INDEXTYPE int
#include <flens/flens.cxx>

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

using namespace std;
using namespace flens;

// typedef mpf_class   T;

typedef double   T;


#define CASE1

int
main()
{
    typedef GeMatrix<FullStorage<T, ColMajor> >   GeMatrix;
    typedef GeMatrix::IndexType                   IndexType;

    typedef DenseVector<Array<T> >                DenseVector;

    const Underscore<IndexType> _;

#   if defined(CASE1)
    const IndexType n = 7;
#   elif defined(CASE2)
    const IndexType n = 4;
#   endif

    GeMatrix            A(n, n), Z(n, n);
    DenseVector         wr(n), wi(n);


#   if defined(CASE1)
    A =  2,   3,  -1,   0,  2,   3,  -1,
        -6,  -5,   0,   2, -6,  -5,   0,
         2,  -5,   6,  -6,  2,  -5,   6,
         2,   3,  -1,   0,  2,   3,  -1,
        -6,  -5,   0,   2, -6,  -5,   0,
         2,  -5,   6,  -6,  2,  -5,   6,
         4,   6,   2,  -3,  4,   6,   2;
#   elif defined(CASE2)
     A =  2,   3,  -1,   0,
         -6,  -5,   0,   2,
          2,  -5,   6,  -6,
          2,   3,  -1,   0;
#   endif

    cerr << "A = " << A << endl;

    DenseVector tau(n-1);
    DenseVector work;

    lapack::hrd(IndexType(1), n, A, tau, work);

    cerr << "A = " << A << endl;
    cerr << "lange(OneNorm, A, work) = "
         << lapack::lange(lapack::OneNorm, A, work)
         << endl;
}
