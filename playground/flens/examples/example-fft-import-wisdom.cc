#include <iostream>
#define USE_PLAYGROUND
#define WITH_FFTW

///
/// Set flag to import wisdom and specify file
/// Irrelevant for Intel MKL (does not support wisdom)
///
#define FFTW_WISDOM_IMPORT
#define FFTW_WISDOM_FILENAME "wisdom.dat"
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

    const IndexType n = 131072;
    Vector x(n), y(n), z(n);

    ///
    /// Fill in random values
    ///
    fillRandom(x);

    ///
    /// Test Fourier Transform with wisdom
    /// Get time via "time ./example-fft-export-wisdom"
    ///

    dft_forward(x, y);

    return 0;
}
