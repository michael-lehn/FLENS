#ifdef USE_SIMPLE_ET
#include <flens/examples/simple_et.h>
#endif

#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 15;
    DenseVector<Array<double> >   x1(n), x2(n), x3(n), y(n);

    x1 = 1;
    x2 = 2;
    x3 = 3;

    y += 2.0*x1 + 4.0*x2 + 3.0*x3;
    cout << "y = " << y << endl;

    y = 2.0*x1 + 4.0*x2 + 3.0*x3;
    cout << "y = " << y << endl;
}
