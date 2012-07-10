
#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>             Complex;
    const Complex                       I(0,1);

    GeMatrix<FullStorage<Complex> >     A(4, 4);
    DenseVector<Array<Complex> >        b(4);
    DenseVector<Array<int> >            piv(4);

///
/// Setup matrix $A$.  Note that `2+I` or `2*I` would not work.  That's because
/// the STL does not define the operation `int+complex<double>` or
/// `int*complex<double>`.
///
    A =   2,    3,  -I,    0,
         -6,   -5,   0, 2.+I,
       2.*I,   -5,   6,   -6,
          4, 6.*I,   2,   -3;

    b = 20,
       -33.*I,
       -43,
        49;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the $LU$ factoriation of $A$.
///
    int info = lapack::trf(A, piv);

///
/// If `info` is not zero then the matrix was singular.  More precise, the
/// upper triangular matrix $U$ computed by `trf` is exactly singular.
///
    if (info!=0) {
        cerr << "A is singular" << endl;
        return info;
    }

///
///  Compute $A^{-1}$ from the $LU$ factorization of $A$.
///
    lapack::tri(A, piv);
    cout << "inv(A) = " << A << endl;

///
/// Compute $x = A^{-1}b$.
///
    DenseVector<Array<Complex> >        x;
    x = A*b;

    cout << "x = inv(A)*b = " << x << endl;

    return 0;
}

