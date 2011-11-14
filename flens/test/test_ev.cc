//
// compile with:
/*
    clang++ -std=c++0x -DWITH_REFBLAS -c -DCHECK_CXXLAPACK test_ev.cc  -I../../
    gfortran test_ev.o -lstdc++ ilaenv.o iparmq.o drivers.o ../lapack/interface/cblas_REF.a ../lapack/interface/blas_REF.a ../lapack/interface/lapack_FLENS.a
*/
//

#include <iostream>
#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


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
#   elif defined(CASE3)
    const IndexType n = 4;
#   endif

    GeMatrix            A(n, n), VL(n, n), VR(n, n);
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
#   elif defined(CASE3)
     A =  1,   9,   1,   2,
          1,   2,   1,   2,
          0,   1,   3,   2,
          1,   0,   1,   4;
#   endif

    cerr << "A = " << A << endl;

    DenseVector work(5*n);

    lapack::ev(false, true, A, wr, wi, VL, VR, work);

    cerr << "A = " << A << endl;
    cerr << "wr = " << wr << endl;
    cerr << "wi = " << wi << endl;
    cerr << "VL = " << VL << endl;
    cerr << "VR = " << VR << endl;
}
