#include <complex>
#include <iostream>

#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
///
/// Some typedefs for real and complex valued matrices.
///
   typedef GeMatrix<FullStorage<complex<double> > >    ComplexGeMatrix;
   typedef GeMatrix<FullStorage<double> >              RealGeMatrix;

    const int m = 5;
    const int n = 4;
    ComplexGeMatrix  Z(m,n);

///
/// See the next page for details on `IndexVariable` and `Complex`.
///
    ComplexGeMatrix::IndexVariable  i,j;
    Z(i,j) = Complex(i+j,j-i);

///
/// You can initialize a real valued matrix through the constructor.
///
    RealGeMatrix X = real(Z);

///
///  Or through an assignment.
///
    RealGeMatrix  Y;
    Y = imag(Z);

    cout << "Z = " << Z << endl;
    cout << "real(Z) = " << X << endl;
    cout << "imag(Z) = " << Y << endl;

///
/// You also can overwrite the real- or imaginary-part of a complex matrix with
/// a real valued matrix.  In this case we swap the real- and imaginary parts
/// of $Z$.  This is possible as $X$ and $Y$ contain copies.
///
    real(Z) = Y;
    imag(Z) = X;

    cout << "Z = " << Z << endl;
}
