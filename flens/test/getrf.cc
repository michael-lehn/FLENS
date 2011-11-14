//
// compile with:
//
// clang++ -std=c++0x -Wall -I../../ getrf.cc
//
//

#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define INTEGER         int
#define DOUBLE          double
#define LOGICAL         int
#define FLOAT           float
#define DOUBLE_COMPLEX  double
#define UNKNOWN         char*

//------------------------------------------------------------------------------
#include <iostream>
#include <flens/flens.cxx>


using namespace std;
using namespace flens;

// typedef mpf_class   T;


typedef double   T;

#define TEST_CASE_3

int
main()
{
    typedef GeMatrix<FullStorage<T> >           GeMatrix;
    typedef GeMatrix::IndexType                 IndexType;

    typedef DenseVector<Array<IndexType> >      IndexVector;

    const Underscore<IndexType> _;

    const IndexType m = 2,
                    n = 3;

    GeMatrix            A(m, n);
    IndexVector         iPiv(m);

    A =   -0.689381,                    0,                    0,
           0.498709,            -0.725288,                    0;

    cerr << "A = " << A << endl;

    lapack::trf(A, iPiv);

    // lapack::sv(Ab(_,_(1,m)), iPiv, Ab(_,_(m+1,n)));
    //lapack::trf(Ab(_,_(1,m)), iPiv);
    //lapack::trs(transA, Ab(_,_(1,m)), iPiv, Ab(_,_(m+1,n)));

    cerr << "-> A = " << A << endl;
    cerr << "-> iPiv = " << iPiv << endl;
}
