#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;


typedef GeMatrix<FullStorage<double> >  DGeMatrix;
typedef DenseVector<Array<double> >     DDenseVector;


// Each function call causes a malloc for a 4x4 matrix and a free at the end of
// scope.
void
foo()
{
    DGeMatrix A(4,4);
    A =  1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12,
        13, 14, 15, 16;

    cout << "A = " << A << endl;
}

int
main()
{
    // Always have in mind that the function is an essential inner core routine
    // doing only little computation but gets called many times.
    for (int i=0; i<100; ++i) {
        foo();
    }
}
