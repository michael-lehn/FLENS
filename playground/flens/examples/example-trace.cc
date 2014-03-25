#include <iostream>
#define USE_PLAYGROUND
#include <flens/flens.cxx>

using namespace std;
using namespace flens;
using namespace dft;

typedef complex<double>   T;

int
main(int argc, char* argv[])
{
    ///
    ///  Define some convenient typedefs for the vector types
    ///
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef Matrix::IndexType                   IndexType;

    const IndexType n = 4;
    Matrix A(n, n);

    /// 
    /// Fill in random values
    ///
    fillRandom(A);
    cout << " A = " << A << endl;
    
    ///
    /// Calculate trace
    ///
    cout << lapack::extensions::trace(A) << endl;
        
    return 0;
}
