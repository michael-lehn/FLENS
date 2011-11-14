//
// compile with:
/*
    clang++ -std=c++0x -DWITH_REFBLAS -c -DCHECK_CXXLAPACK test_hd2.cc  -I../../
    gfortran test_hd2.o -lstdc++ ilaenv.o iparmq.o drivers.o ../lapack/interface/cblas_REF.a ../lapack/interface/blas_REF.a ../lapack/interface/lapack_FLENS.a
*/
//

#include <iostream>

#include <flens/lapack/interface/include/config.h>

#include <flens/flens.cxx>

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

using namespace std;
using namespace flens;

// typedef mpf_class   T;

typedef double   T;

#define CASE4

int
main()
{
    typedef GeMatrix<FullStorage<T, ColMajor> >     GeMatrix;
    typedef GeMatrix::IndexType                     IndexType;

    typedef DenseVector<Array<T> >                  DenseVector;

    const Underscore<IndexType> _;

    const IndexType n = 4;

    GeMatrix            A(n, n), Z(n, n);
    DenseVector         wr(n), wi(n);

#   if defined(CASE1)

    A =  3,  -1,   0,   0,
        -1,   3,  -1,   0,
         0,  -1,   3,  -1,
         0,   0,  -1,   3;

#   elif defined(CASE2)

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

#   elif defined(CASE3)

    A =  2,  -6,   2,   4,
         3,  -5,  -5,   6,
        -1,   0,   6,   2,
         0,   2,  -6,  -3;

#   elif defined(CASE4)

     A =  1,   9,   1,   2,
          1,   2,   1,   2,
          0,   1,   3,   2,
          1,   0,   1,   4;

#endif

    cerr << "A = " << A << endl;

    DenseVector tau(n-1);
    DenseVector work(n);

    lapack::hrd(IndexType(1), n, A, tau, work);

    cerr << "A = " << A << endl;
    /*
    cerr << "lange(OneNorm, A, work) = "
         << lapack::lange(lapack::OneNorm, A, work)
         << endl;
    */
/*
    IndexType info = lapack::hseqr(lapack::HSEQR::Eigenvalues, lapack::HSEQR::Init,
                                   IndexType(1), n, A, wr, wi, Z, work);

    cerr << "-> info =  " << info << endl;
    cerr << "-> A =  " << A << endl;
    cerr << "-> wr = " << wr << endl;
    cerr << "-> wi = " << wi << endl;
    cerr << "-> Z =  " << Z << endl;
*/
}
