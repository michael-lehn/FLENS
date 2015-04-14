#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >      DGeMatrix;
    typedef DGeMatrix::ConstView                          DGeMatrixConstView;

    const Underscore<DGeMatrix::IndexType>  _;

    DGeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

    DGeMatrixConstView  B = A(_(2,3),_);

    cout << "B(1,1) = " << B(1,1) << endl;
}
