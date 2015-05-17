#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;


typedef GeMatrix<FullStorage<double> >  DGeMatrix;
typedef DenseVector<Array<double> >     DDenseVector;

int
main()
{
    DGeMatrix     A(3,3);

    fillRandom(A);

    cout << "A = " << A << endl;

    cout << "A.strideRow()          = " << A.strideRow() << endl;
    cout << "A.strideCol()          = " << A.strideCol() << endl << endl;

    cout << "A.data()               = " << A.data() << endl;
    cout << "&A(1,1)                = " << &A(1,1) << endl << endl;

    cout << "A.data()+A.strideRow() = " << A.data()+A.strideRow() << endl;
    cout << "&A(2,1)                = " << &A(2,1) << endl << endl;

    cout << "A.data()+A.strideCol() = " << A.data()+A.strideCol() << endl;
    cout << "&A(1,2)                = " << &A(1,2) << endl << endl;
}
