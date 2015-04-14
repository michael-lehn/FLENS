#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >  DGeMatrix;
    typedef DenseVector<Array<double> >     DDenseVector;

    Underscore<DDenseVector::IndexType>     _;

    ///
    ///  Allocate and initialize some vector ...
    ///
    DDenseVector x(5);
    x = 1, 2, 3, 4, 5;

    ///
    ///  ... and matrix.
    ///
    DGeMatrix A(4, 5);
    A = 11, 12, 13, 14, 15,
        21, 22, 23, 24, 25,
        31, 32, 33, 34, 35,
        41, 42, 43, 44, 45;

    ///
    ///  Print the vector and matrix
    ///
    cout << "x = " << x << endl;
    cout << "A = " << A << endl;
    cout << endl;

    ///
    ///  Create a vector view that references every second element of `x`.
    ///
    DDenseVector::View  y = x(_(1,2,5));

    cout << "y = x(1:2:5) = " << y << endl;
    cout << endl;


    ///
    ///  Create a vector view that references the second row of `A`
    ///
    DDenseVector::View row = A(2,_);

    cout << "row = A(2,:) = " << row << endl;
    cout << endl;

    ///
    ///  Create a vector view that references part of second row of `A`
    ///
    DDenseVector::View rowPart = A(2,_(2,4));

    cout << "rowPart = A(2,2:4) = " << rowPart << endl;
    cout << endl;


    ///
    ///  Create a vector view that references the third column of `A`
    ///
    DDenseVector::View col = A(_,3);

    cout << "col = A(:,3) = " << col << endl;
    cout << endl;

    ///
    ///  Create a vector view that references the diagonal of `A`
    ///
    DDenseVector::View d0 = A.diag(0);

    cout << "d0 = A.diag(0) = " << d0 << endl;
    cout << endl;

    ///
    ///  Create a vector view that references the first sub-diagonal of `A`
    ///
    DDenseVector::View d_1 = A.diag(-1);

    cout << "d_1 = A.diag(-1) = " << d_1 << endl;
    cout << endl;

    ///
    ///  Create a vector view that references the first super-diagonal of `A`
    ///
    DDenseVector::View d1 = A.diag(1);

    cout << "d1 = A.diag(1) = " << d1 << endl;
    cout << endl;


    return 0;
}
