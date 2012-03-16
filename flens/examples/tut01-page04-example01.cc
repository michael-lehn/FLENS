#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{

///
/// The storage scheme for `DenseVector` is `Array`, `ArrayView` or
/// `ConstArrayView`.  We setup a typedef for a "regular" dense vector that
/// allocates its own memory.
///
    typedef DenseVector<Array<double> >   DenseVector;
    typedef DenseVector::IndexType        IndexType;

///
/// We define a dense vector of length `4` and initialize it.
///
    DenseVector x(4);
    x = 1, 2, 3, 4;

///
/// We retrieve some information like index range and vector length.  You
/// will see that the default index base is One.
///
    cout << "x.range() = " << x.range() << endl;
    cout << "x.length() = " << x.length() << endl;

///
/// We print the vector.
///
    cout << "x = " << x << endl;

///
/// We iterate through the vector and reset its elements (such that $x_i =
/// i^2$).
///
    for (IndexType i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        x(i) = i*i;
    }

///
/// For selecting vector parts we define a range operator
///
    const Underscore<IndexType> _;

///
/// We create a vector view `y` that references a part of `x`.
///
    DenseVector::View y = x(_(2,4));
    y = 666;

///
/// Here we create a new vector `z` that allocates its own memory and
/// initializes it with a copy of a part of `x`.  So `z` is no vector view.
///
    DenseVector::NoView z = x(_(1,2));
    z = 42;

    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
}
