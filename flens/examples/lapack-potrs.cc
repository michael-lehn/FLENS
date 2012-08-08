#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(3, 3);
    DenseVector<Array<double> >      b(3);

///
/// Setup the raw data.
///
    A = 2, 1, 0,
        0, 2, 1,
        0, 0, 2;

    b = 1,
        2,
        3;

///
/// $S$ is a symmetric matrix view constructed from the upper triangular
/// part of $A$.  Note that $S$ only references data from $A$.
///
    auto S = A.upper().symmetric();

    cout << "S = " << S << endl;
    cout << "b = " << b << endl;

///
/// Computes the Cholesky factorization $S = L L^T$ where $L$ is upper
/// triangular.  Matrix $S$ (i.e. the lower triangular part of $A$) gets
/// overwritten with $L^T$.
///
    int info = lapack::potrf(S);

///
/// If `info` is not zero the factorization could not be computed as the matrix
/// was singular.
///
    if (info!=0) {
        cerr << "S is singular" << endl;
        return info;
    }

///
/// Compute $x = S^{-1}b$.  Vector $b$ gets overwritten with the solution
/// vector $x$.
///
    lapack::potrs(S, b);
    cout << "inv(S)*b = " << b << endl;

    return 0;
}

