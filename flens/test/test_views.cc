//
// compile with:
//
// clang++ -std=c++0x test_views.cc -I../../
//
//

#include <iostream>
#include <flens/flens.cxx>

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

using namespace std;
using namespace flens;

typedef double  FLOAT;
typedef double  CXXFLOAT;
typedef int     INT;


typedef IndexOptions<INT>       I;
typedef std::allocator<FLOAT>   Alloc;

typedef FullStorageView<CXXFLOAT, cxxblas::ColMajor, I, Alloc>  FSV;
typedef FullStorage<CXXFLOAT, cxxblas::ColMajor>                FS;


typedef ArrayView<CXXFLOAT, I, Alloc>                           AV;
typedef ArrayView<INT, I, Alloc>                                IAV;


int
main()
{
    INT     M = 4, N= 3;
    INT     n = M*N;

    FLOAT   x[n];
    INT     iPiv[n];

    for (INT i=0; i<n; ++i) {
        x[i] = i+1;
        iPiv[i] = i;
    }

    DenseVector<AV> y = AV(x, Alloc(), n, INT(1));
    DenseVector<IAV> piv = IAV(iPiv, Alloc(), n, INT(1));
    GeMatrix<FSV>   A = FSV(x, Alloc(), M, N, 2);
    GeMatrix<FS>    B = GeMatrix<FSV>(FSV(x, Alloc(), M, N, 2));


    cerr << "y = " << y << endl;
    cerr << "piv = " << piv << endl;
    cerr << "piv.firstIndex() = " << piv.firstIndex() << endl;
    
    cerr << "A = " << A << endl;
    cerr << "B = " << B << endl;
    y(2) = 5;

    cerr << "x = " << endl;
    for (INT i=0; i<n; ++i) {
        cerr << x[i] << " ";
    }
    cerr << endl;
    cerr << "A = " << A << endl;
    cerr << "B = " << B << endl;
}
