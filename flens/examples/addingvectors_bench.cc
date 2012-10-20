#ifdef USE_SIMPLE_ET
#include <flens/examples/simple_et.h>
#elif defined USE_SIMPLE_ET2
#include <flens/examples/simple_et2.h>
#endif

#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    const int n = 150;
    DenseVector<Array<double> >   x1(n), x2(n), x3(n), y(n);

    x1 = 1;
    x2 = 2;
    x3 = 3;

    for (int count=1; count<=5000000; ++count) {
        y += 2.0*x1 + 4.0*x2 + 3.0*x3;
        y = 2.0*x1 + 4.0*x2 + 3.0*x3;
    }
}
