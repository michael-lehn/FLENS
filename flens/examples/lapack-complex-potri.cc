#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>           Complex;
    const Complex                     I(0,1);

    GeMatrix<FullStorage<Complex> >   A(3, 3);

///
/// Setup the raw data.
///
    A = 2,    0, 0,
        I,    2, 0,
        0, 1.+I, 2;

///
/// $H$ is a symmetric matrix view constructed from the lower triangular
/// part of $A$.  Note that $H$ only references data from $A$.
///
    auto H = A.lower().hermitian();
    cout << "H = " << H << endl;

///
/// Computes the Cholesky factorization $H = L L^T$ where $L$ is lower
/// triangular.  Matrix $H$ (i.e. the lower triangular part of $A$) gets
/// overwritten with $L$.
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
/// Compute $H^{-1}$ from the Cholesky factorization $H = L L^T$.
///
    cout << "H = " << H << endl;
    lapack::potri(H);
    cout << "inv(H) = " << H << endl;

    return 0;
}

