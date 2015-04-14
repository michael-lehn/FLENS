#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double, ColMajor> >      DGeMatrix;
    typedef DGeMatrix::View                               DGeMatrixView;
    typedef DGeMatrix::ConstView                          DGeMatrixConstView;

    const Underscore<DGeMatrix::IndexType>  _;

    DGeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

    auto       B = A(_(2,3),_);       // auto: DGeMatrixView,
                                      // B:    DGeMatrixView

    const auto C = B(_,_(2,4));       // auto: DGeMatrixView,
                                      // C:    const DGeMatrixView

    const auto D = C;                 // auto: DGeMatrixView,
                                      // D:    const DGeMatrixView

    const auto E = C(_,_(1,2));       // auto: DGeMatrixConstView,
                                      // E: const DGeMatrixConstView

    // auto F = C(_,_(1,2));          // auto: DGeMatrixConstView
                                      // F:    DGeMatrixConstView
    // cout << "F(1,1) = " << F(1,1)  // error: static assertion failure as you
    //      << endl;                  //        call the non-const operator.

}
