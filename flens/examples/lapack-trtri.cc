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

///
/// $T$ is an upper triangular matrix referencing the upper triangular part
/// of matrix $A$.
///
    auto T = A.upper();
    cout << "T = " << T << endl;

///
/// Computes the inverse $T^{-1}$ of matrix $T$.
///
    int info = lapack::tri(T);

///
/// If `info` is not zero then matrix $T$ was singular.  To be more precise, in
/// this case one diagonal element was exactly zero.
///
    if (info!=0) {
        cerr << "T is singular" << endl;
        return info;
    }

    cout << "inv(T) = " << T << endl;

    return 0;
}

