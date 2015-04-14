#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;


typedef GeMatrix<FullStorage<double> >  DGeMatrix;
typedef DenseVector<Array<double> >     DDenseVector;


void
foo()
{
    double buffer_[4*4] = { 1,  2,  3,  4,
                            5,  6,  7,  8,
                            9, 10, 11, 12,
                           13, 14, 15, 16};

    DGeMatrix::View    A = DGeMatrix::EngineView(4, 4, buffer_, 1, 4);
    DGeMatrix::View    B = DGeMatrix::EngineView(4, 4, buffer_, 4, 1);
    DGeMatrix::View    C = DGeMatrix::EngineView(2, 4, buffer_, 4, 1);
    DDenseVector::View x = DDenseVector::EngineView(4*4, buffer_);
    DDenseVector::View y = DDenseVector::EngineView(2*4, buffer_, 2);

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
}

int
main()
{
    for (int i=0; i<100; ++i) {
        foo();
    }
}
