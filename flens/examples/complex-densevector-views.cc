#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>        Complex;
    DenseVector<Array<Complex> >   z(5);

    z = Complex(1,2), Complex(3,4), Complex(5,6), Complex(7,8), Complex(9,10);

    Underscore<int> _;

///
/// $x$ and $y$ are of type `DenseVector<ArrayView<double> >` but `auto` is
/// shorter.
///
    auto x = real(z);
    auto y = imag(z);

    cout << "z = " << z << endl;
    cout << "real(z) = " << x << endl;
    cout << "imag(z) = " << y << endl;

///
/// Let's modify the real and imaginary parts of $z$ through the vector views.
/// Note, $x$ and $y$ are of type `DenseVector`.  So you just can use them as
/// any other `DenseVector`.
///
    x = 2, 9, 4, 7, 6;
    y(_(1,2,5)) = 666, -666, 666;

///
/// Now we check that $z$ actually was modified:
///
    cout << "z = " << z << endl;

///
/// As $x$ and $y$ are dense vectors you can apply any BLAS function
///
    blas::swap(x,y);
    cout << "z = " << z << endl;
}
