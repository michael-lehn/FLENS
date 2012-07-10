#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    GeMatrix<FullStorage<double> >   A(3,3);
    DenseVector<Array<int> >         piv(3);

///
/// Setup the toy example.
///
    A = 1, 2, 3,
        4, 5, 6,
        7, 8, 1;

    cout << "A = " << A << endl;

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

    return 0;
}

