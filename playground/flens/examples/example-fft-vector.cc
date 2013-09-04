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
    typedef DenseVector<Array<T> >              Vector;
    typedef Vector::IndexType                   IndexType;

    const IndexType n = 4;
    Vector x(n), y(n), z(n);

    /// 
    /// Fill in random values
    ///
    fillRandom(x);
        
    ///
    /// Calculate forward Fourier transform
    ///
    dft_forward(x, y);
    
    ///
    /// Calculate normalized backward Fourier transform
    ///
    dft_backward_normalized(y, z);
    
    ///
    /// Check results
    ///
    cout << " x = " << x << endl;
    cout << " dft(x) = " << y << endl;
    cout << " dft^{-1}(dft(x)) = " << z << endl;

    return 0;
}
