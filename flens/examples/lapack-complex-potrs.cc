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
    A = 2, I,    0,
        0, 2, 1.+I,
        0, 0,    3;

    b = 1,
        2,
        3;

///
/// $H$ is a symmetric matrix view constructed from the upper triangular
/// part of $A$.  Note that $H$ only references data from $A$.
///
    auto H = A.upper().hermitian();

    cout << "H = " << H << endl;
    cout << "b = " << b << endl;

///
/// Computes the Cholesky factorization $H = L L^T$ where $L$ is upper
/// triangular.  Matrix $H$ (i.e. the lower triangular part of $A$) gets
/// overwritten with $L^T$.
///
    int info = lapack::potrf(H);

///
/// If `info` is not zero the factorization could not be computed as the matrix
/// was singular.
///
    if (info!=0) {
        cerr << "H is singular" << endl;
        return info;
    }

///
/// Compute $x = H^{-1}b$.  Vector $b$ gets overwritten with the solution
/// vector $x$.
///
    lapack::potrs(H, b);
    cout << "inv(H)*b = " << b << endl;

    return 0;
}

