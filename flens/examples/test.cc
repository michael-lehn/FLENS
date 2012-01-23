#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef DenseVector<Array<T> >      DenseVector;

    DenseVector    a(n), b(n), c(n), d(n);
    IndexVector    piv(n);

    a = 1;
    b = 2;
    c = 3;
    d = 4;

    d = a + b + c;

    cerr << "d = " << d << endl;
}
