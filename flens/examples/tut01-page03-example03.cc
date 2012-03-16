#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

///
/// This function receives a general matrix and `MA` can be of type
/// `FullStorage`, `FullStorageView` or `ConstFullStorageView`. Note that the
/// function receives the matrix as a const reference.
template <typename MA>
void
dummy(const GeMatrix<MA> &A)
{

///
/// We setup a range operator for this matrix type ...
///
    typedef typename GeMatrix<MA>::IndexType  IndexType;
    const Underscore<IndexType>  _;

///
/// ... and get some numbers we will use to split the matrix vertically.
///
    const IndexType n = A.numCols();
    const IndexType k = n/2;

///
/// As `A` is in this function in const context the created matrix views will
/// always be of type `GeMatrix<ConstFullStorageView<..> >`. And this is also
/// the type represented by `auto`. So `A1` and `A2` will be of type
/// `const GeMatrix<ConstFullStorageView<..> >`.
///
    const auto A1 = A(_,_(1,k));
    const auto A2 = A(_,_(k+1,n));

    cout << "A1 = " << A1 << endl;
    cout << "A2 = " << A2 << endl;
}

///
/// In the `main` function we again start with a typedef for a `GeMatrix` with
/// elements of type `double` that are stored column-wise in memory. Also, we
/// define a suitable range operator.
///
int
main()
{
    typedef FullStorage<double, ColMajor>  FS;
    typedef GeMatrix<FS>                   GeMatrix;
    const Underscore<GeMatrix::IndexType>  _;

///
/// Then we setup some matrix.
///
    GeMatrix  A(3,4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12;

    cout << "A = " << A << endl;

///
/// Now we call dummy with different view types.
///
/// Call it with matrix of type `GeMatrix<FullStorage<..> >`
///
    dummy(A);

///
/// Call it with matrix of type `GeMatrix<FullStorageView<..> >`
///
    dummy(A(_(1,3),_));

///
/// Call it again with matrix of type `GeMatrix<FullStorageView<..> >`
///
    auto B = A(_(1,3),_);
    dummy(B);

///
/// Call it with matrix of type `GeMatrix<ConstFullStorageView<..> >`. Here
/// `C` is of type `const GeMatrix<FullStorageView>`, i.e. a (non-const) matrix
/// view in const context.
///
    const auto C = A(_(1,3),_);
    dummy(C(_,_(1,4)));
}
