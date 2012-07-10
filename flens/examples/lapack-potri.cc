#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(3, 3);

///
/// Setup the raw data.
///
    A = 2, 0, 0,
        1, 2, 0,
        0, 1, 2;

///
/// $S$ is a symmetric matrix view constructed from the lower triangular
/// part of $A$.  Note that $S$ only references data from $A$.
///
    auto S = A.lower().symmetric();
    cout << "S = " << S << endl;

///
/// Computes the Cholesky factorization $S = L L^T$ where $L$ is lower
/// triangular.  Matrix $S$ (i.e. the lower triangular part of $A$) gets
/// overwritten with $L$.
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
/// Compute $S^{-1}$ from the Cholesky factorization $S = L L^T$.
///
    lapack::potri(S);
    cout << "inv(S) = " << S << endl;

    return 0;
}

