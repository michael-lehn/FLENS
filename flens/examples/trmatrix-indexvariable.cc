#include <complex>
#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >             RealGeMatrix;
    typedef TrMatrix<FullStorage<double> >             RealTrMatrix;
    typedef SyMatrix<FullStorage<double> >             RealSyMatrix;
    typedef HeMatrix<FullStorage<complex<double> > >   ComplexHeMatrix;
    typedef GeMatrix<FullStorage<complex<double> > >   ComplexGeMatrix;

///
/// We allocate a triangular matrix, define index variables and use them.
/// As the matrix is upper triangular $i$ and $j$ will only iterate over the
/// upper triangular part.
///
    RealTrMatrix   A(3, Upper);
    RealTrMatrix::IndexVariable i, j;
    A(i,j) = 10*i+j;

    cout << "A = " << A << endl;


///
/// For symmetric and hermitian matrices ...
///
    RealSyMatrix      B(3, Upper);
    ComplexHeMatrix   C(3, Lower);

///
/// ... you simply initialize their dominant triangular part
///
    B.triangular()(i,j) = 2+i-j;
    C.triangular()(i,j) = Complex(i+j,j);

    cout << "B = " << B << endl;
    cout << "C = " << C << endl;

///
/// Finally we initialize different triangular parts of a general complex
/// valued matrix
///
    ComplexGeMatrix   D(3,4);
    D.strictLower()(i,j) = Complex(i,0);
    D.diag(0)(i)         = Complex(666,i);
    D.strictUpper()(i,j) = Complex(i,666);

    cout << "D = " << D << endl;
}
