#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef complex<double>                     ZComplex;
    typedef GeMatrix<FullStorage<ZComplex> >    ZGeMatrix;

    ZComplex   I(0,1);

    ZGeMatrix  A(2,3);

    A = 2.*I, 3.+I,    1.,
        4.,   2.+2.*I, 4.;

    cout << "A = " << A << endl;

    // in-place conjugate (via overloaded operators)
    A = conjugate(A);
    cout << "A = " << A << endl;

    // in-place conjugate (via explicit copy call)
    blas::copy(Conj, A, A);
    cout << "A = " << A << endl;

    // in-place conjugate (via explicit call)  (new)
    blas::conj(A);
    cout << "A = " << A << endl;

    // elementwise variant
    ZGeMatrix::IndexVariable  k, l;

    A(k,l) = Real(A(k,l)) - Imag(A(k,l)) * I;
    cout << "A = " << A << endl;
}
