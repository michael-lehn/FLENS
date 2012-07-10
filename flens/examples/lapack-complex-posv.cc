#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>             Complex;
    const Complex                       I(0,1);

    GeMatrix<FullStorage<Complex> >     A(3, 3);
    DenseVector<Array<Complex> >        b(3);

///
/// Setup the raw data.
///
    A = 2,    0, 0,
        I,    2, 0,
        0, 1.+I, 3;

    b =    1,
        2.+I,
           3;

///
/// $H$ is a hermitian matrix view constructed from the lower triangular
/// part of $A$.  Note that $H$ only references data from $A$.
///
    auto H = A.lower().hermitian();

    cout << "H = " << H << endl;
    cout << "b = " << b << endl;

///
/// - Computes the Cholesky factorization $H = L L^T$ where $L$ is lower
///   triangular.  Matrix $H$ (i.e. the lower triangular part of $A$) gets
///   overwritten with $L$.
/// - Solves $Lu=b$ and then $L^T x = u$.  Vector $b$ gets overwritten with $x$.
///
    int info = lapack::posv(H, b);

///
/// If `info` is not zero the factorization could not be computed as the matrix
/// was singular.
///
    if (info!=0) {
        cerr << "H is singular" << endl;
        return info;
    }

    cout << "inv(H)*b = " << b << endl;

    return 0;
}

