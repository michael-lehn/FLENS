#include <complex>
#include <iostream>

#include <flens/examples/complexhack.h>
#include <flens/examples/complexhack.tcc>
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
    auto x = real(z);
    auto y = imag(z);

    cout << "z = " << z << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    //  let's modify the real part of z

    x(_(2,4)) = -1, -5, -17;

    //  As x is a view we changed the original vector
    cout << "We changed x:" << endl;
    cout << "z = " << z << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    //  imaginary numbers are evil!
    y(_(1, 2, 5)) = 666, -666, 666;

    cout << "We changed z:" << endl;
    cout << "z = " << z << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    //  Setup some complex valued matrix Z
    const int m = 5;
    const int n = 4;
    GeMatrix<FullStorage<Complex> >  Z(m,n);

    for (int j=1; j<=n; ++j) {
        for (int i=1; i<=m; ++i) {
            Z(i,j) = Complex(i,j);
        }
    }

    //
    // Now we show how you can use real(Z) and imag(Z) to copy real and
    // imaginary part to real valued matrices.  Note that in the following
    // matrices X and Y must not be views.
    //

    {
        //  You can initialize a real valued matrix through the constructor
        GeMatrix<FullStorage<double> >  X = real(Z);

        //  Or you can assign it.
        GeMatrix<FullStorage<double> >  Y;
        Y = imag(Z);

        cout << "Z = " << Z << endl;
        cout << "X = " << X << endl;
        cout << "Y = " << Y << endl;

        //
        //  You also can copy things back.  As X and Y contain copies we can
        //  swap real and imaginary parts.
        //
        real(Z) = Y;
        imag(Z) = X;

        cout << "Z = " << Z << endl;
    }

    {
        // You also can transpose a copy
        GeMatrix<FullStorage<Complex> >  W;
        GeMatrix<FullStorage<double> >   X, Y;

        X = 3.0*transpose(real(Z));
        Y = imag(Z);

        cout << "X = " << X << endl;
        cout << "Y = " << Y << endl;

        imag(W) = X;
        real(W) = transpose(Y);
        cout << "W = " << W << endl;
    }

    {
        // Let's check const-correctness
        const GeMatrix<FullStorage<Complex> >  V = Z;
        GeMatrix<FullStorage<double> >   X, Y;

        cout << "V = " << V << endl;

        X = real(V);
        Y = imag(V);

        cout << "X = " << X << endl;
        cout << "Y = " << Y << endl;

        // if you uncomment this compilation should fail:
        //real(V) = X;
        //imag(V) = Y;
    }
}
