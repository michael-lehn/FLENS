#include <cxxstd/iostream.h>

#define USE_PLAYGROUND
#define WITH_FFTW

///
/// Possbile options are:
///   WITH_FFTW_FLOAT, WITH_FFTW_DOUBLE (default),
///   WITH_FFTW_LONGDOUBLE, WITH_FFTW_QUAD
/// Link it against
/// -lfftwf        , -lfftw                    , -lfftwl             , -lfftwq
///
#define WITH_FFTW_FLOAT

#include <flens/flens.cxx>

using namespace std;
using namespace flens;
using namespace dft;

typedef complex<float>   T;

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
