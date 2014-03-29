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

    const IndexType n = 4,
                    m = 8;
    Matrix A(n, m), B(n, m), C(n, m);

    ///
    /// Fill in random values
    ///
    fillRandom(A);

    ///
    /// Calculate forward Fourier transform
    /// of every column of A.
    /// Use dft_row_forward for rows.
    ///
    dft_col_forward(A, B);

    ///
    /// Calculate normalized backward Fourier transform
    /// of every column of A
    /// Use dft_row_backward_normalized for rows
    ///
    dft_col_backward_normalized(B, C);

    ///
    /// Check results
    ///
    cout << " A = " << A << endl;
    cout << " dft(A) = " << B << endl;
    cout << " dft^{-1}(dft(A)) = " << C << endl;

    return 0;
}
