#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
///
/// We define typedefs for dense vectors and general matrices that allocate
/// their own memory.
///
    typedef GeMatrix<FullStorage<double, ColMajor> >    GeMatrix;
    typedef DenseVector<Array<double> >                 DenseVector;
    typedef DenseVector::IndexType                      IndexType;

///
/// For selecting vector parts we define a range operator
///
    const Underscore<IndexType> _;

///
/// We setup some matrix.
///
    GeMatrix  A(4, 5);
    A =  1,  2,  3,  4,  5,
         6,  7,  8,  9, 10,
        11, 12, 13, 14, 15,
        16, 17, 18, 19, 20;

///
/// We create a vector view `y` that references the second row of `A`.  Then
/// we reset all elements of `y` to `666`.  This will also change elements in
/// `A`.
///
    DenseVector::View y = A(2,_);
    y = 666;

///
/// We create a "regular" vector `z` and initialize it with the third row
/// of `A`.  Just to show that `z` is not a view we then reset all elements
/// of `z` to `42`.  You will see that this will not change elements in `A`.
///
    DenseVector::NoView z = A(3,_);
    z = 42;

    cerr << "A = " << A << endl;
    cerr << "y = " << y << endl;
    cerr << "z = " << z << endl;

///
/// Of course we also can use the `auto` type instead of writing out the
/// `DenseVector::View` type.
///
    auto x = A(_(2,4),2);
    x = 999;

///
/// You also can reference the `k`-th off-diagonal of `A`.  For $k=0$ you
/// get the main diagonal, for $k>0$ the `k`-th super-diagonal and for
/// $k<0$ the corresponding sub-diagonal.
///
    auto d_1 = A.diag(-1);
    auto d0  = A.diag(0);
    auto d1  = A.diag(1);
    auto d2  = A.diag(2);

    cerr << "d_1 = " << d_1 << endl;
    cerr << "d0  = " << d0 << endl;
    cerr << "d1  = " << d1 << endl;
    cerr << "d2  = " << d2 << endl;

///
/// Next we set all elements on the diagonal to `-1` and all elements on the
/// first super-diagonal to `-100`, ..., `-400`:
///
    d0 = -1;
    d1 = -100, -200, -300, -400;

    cerr << "A = " << A << endl;
    cerr << "x = " << x << endl;

///
/// We finally print the strides of the vectors.  Compare strides of vectors
/// referencing rows, columns and diagonals.
///
    cerr << "A.leadingDimension() = " << A.leadingDimension() << endl;
    cerr << "x.stride() = " << x.stride() << endl;
    cerr << "y.stride() = " << y.stride() << endl;
    cerr << "z.stride() = " << z.stride() << endl;
    cerr << "d_1.stride() = " << d_1.stride() << endl;
    cerr << "d0.stride()  = " << d0.stride() << endl;
    cerr << "d1.stride()  = " << d1.stride() << endl;
    cerr << "d2.stride()  = " << d2.stride() << endl;
}
