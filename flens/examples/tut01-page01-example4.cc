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

    ///
    ///  We want to create a const-view.  However, we do not declare `B` as
    ///  const.  Ideally the instantiation here would already cause a compile
    ///  time error.  You see next, why.
    ///
    DGeMatrixConstView  B = A(_(2,3),_);

    ///
    ///  `GeMatrix` implements two function operators for element access.  For
    ///  a const-view merely the read-only access must be allowed.  But as `B`
    ///  is not declared as const the read-write access gets called.
    ///
    cout << "B(1,1) = " << B(1,1) << endl;

}
