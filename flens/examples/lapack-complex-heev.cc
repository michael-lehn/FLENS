#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;


int
main()
{
    typedef complex<double>                             ZDouble;
    typedef GeMatrix<FullStorage<ZDouble, ColMajor> >   ZGeMatrix;
    typedef DenseVector<Array<ZDouble> >                ZDenseVector;
    typedef DenseVector<Array<double> >                 DDenseVector;

    const int n = 4;

    ZGeMatrix      A(n, n);
    DDenseVector   w(n);


    A.upper() = ZDouble(1,0), ZDouble(1,1), ZDouble(2,1), ZDouble(5,2),
                              ZDouble(2,0), ZDouble(1,4), ZDouble(2,7),
                                            ZDouble(3,0), ZDouble(1,4),
                                                          ZDouble(4,0);

    cerr << "A.upper().hermitian() = "
         << A.upper().hermitian() << endl;

    ZDenseVector    work;
    DDenseVector    rWork;
    lapack::ev(true, A.upper().hermitian(), w, work, rWork);

    cerr << "A = " << A << endl;
    cerr << "w = " << w << endl;
}
