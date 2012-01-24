#include <iostream>
#include <flens/debug/aux/aux.h>
#include <flens/flens.cxx>
#include <flens/debug/aux/aux.tcc>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef DenseVector<Array<T> >      DenseVector;

    int n = 3;
    DenseVector    a(n), b(n), c(n), d(n), e(n);

    a = 1;
    b = 2;
    c = 3;
    d = 4;
    e = 5;

    verbose::ClosureLog::start("log");

    e = 3*2*a - b + c + d;
    cerr << "e = " << e << endl;

    e = ((a - b) - c) - d;
    cerr << "e = " << e << endl;

    e = a - (b - (c + d));
    cerr << "e = " << e << endl;

    e = a + (b + (c + d));
    cerr << "e = " << e << endl;
}
