#include <complex>
#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >             RealGeMatrix;

///
/// We allocate a matrix and define index variables.
///
    RealGeMatrix   A(3,4);
    RealGeMatrix::IndexVariable i, j;

///
/// In the following statement $i$ and $j$ will automatically iterate depending
/// on the matrix size of $A$.  Obviously $i$ and $j$ are not some plain integer
/// values.  But you can do some simple computations with them.
///
    A(i,j) = 10*i + j;


///
/// Another example which is using an `Abs` function that can operate on
/// such index variables.
///
    RealGeMatrix   B(4,3);
    B(i,j) = Abs(i+j)*(i-j);

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

///
/// Now we also can do some simple element-wise operation.  Note that the
/// range of $i$ and $j$ is determined by the matrix size of $C$
///
    RealGeMatrix   C(3,4);
    C(i,j) = A(i,j) * Sin(M_PI*B(j,i)/4);

    cout << "C = " << C << endl;
}
