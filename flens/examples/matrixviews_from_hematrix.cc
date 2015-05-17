#include <flens/flens.cxx>
#include <complex>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    typedef std::complex<double>      zcomplex;
    GeMatrix<FullStorage<zcomplex> >  A(3,3);

    A =  zcomplex(1,1), zcomplex(1,2), zcomplex(1,3),
         zcomplex(2,1), zcomplex(2,2), zcomplex(2,3),
         zcomplex(3,1), zcomplex(3,2), zcomplex(3,3);

    auto H = A.upper().hermitian();

    cout << "H = " << H << endl;
    cout << "H.triangular() = " << H.triangular() << endl;
    cout << "H.general() = " << H.general() << endl;
    cout << "H.symmetric() = " << H.symmetric() << endl;
}
