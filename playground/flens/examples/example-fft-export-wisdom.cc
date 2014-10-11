#include <cxxstd/iostream.h>

#define USE_PLAYGROUND
#define WITH_FFTW

///
/// Possible options for FFTW_PLANNER_FLAG
/// FFTW_ESTIMATE (default), FFTW_MEASURE, FFTW_PATIENT, ... (see FFTW3 manual)
/// These options are irrelevant if using the Intel MKL
///

#define FFTW_PLANNER_FLAG FFTW_PATIENT

///
/// Set flag to export wisdom and specify file
/// Irrelevant for Intel MKL (does not support wisdom)
///
#define FFTW_WISDOM_EXPORT
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
    /// Create wisdom of a Fourier transform
    ///
    dft_forward(x, y);

    return 0;
}
