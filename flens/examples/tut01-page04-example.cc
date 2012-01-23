#include <flens/flens.cxx>
#include <iostream>

using namespace flens;
using namespace std;

int
main()
{
    ///
    ///  So this is our plain C-array
    ///
    double data[4*4] = { 1,  2,  3,  4,
                         5,  6,  7,  8,
                         9, 10, 11, 12,
                        13, 14, 15, 16};

    ///
    ///  We define one typedef for a storage scheme that only references data
    ///  and another typedef for a corresponding general matrix type. The
    ///  trick is that you can construct `GeMatrix` objects from full storage
    ///  objects.  So what we later do is building a storage object that
    ///  just keeps internally a pointer to the C-array.
    ///
    typedef FullStorageView<double, ColMajor>  FSView;
    typedef GeMatrix<FSView>                   GeMatrixView;

    ///
    ///  We finally create a matrix view that references the C-array as follows:
    ///    - We wrap the C-array in an anonymous full storage view object
    ///    - and construct with it a general matrix view.
    ///  The syntax for the full storage view is
    ///  `FSView(numRows, numCols, data-pointer, leading dimension)`.
    ///
    GeMatrixView  A = FSView(4, 4, data, 4);

    cout << "A = " << A << endl;

    ///
    ///  With a little pointer arithmetic we also can create views for a
    ///  sub-matrix directly (note that `data+5`points to the element of
    ///  `A(2,3)`.
    ///
    GeMatrixView  B = FSView(3, 2, data+5, 4);

    cout << "B = " << B << endl;

    ///
    ///  Note that in the following the data of the C-array gets copied into
    ///  matrix `M`.  This is because `GeMatrixNoView` has a non-view storage
    ///  scheme.
    ///
    typedef GeMatrix<FullStorage<double, ColMajor> >  GeMatrixNoView;
    GeMatrixNoView  M = GeMatrixView(FSView(4, 4, data, 4));

    ///
    ///  Note that you only can construct a `GeMatrix` from a storage object
    ///  if the type of the storage object is identical with `GeMatrix::Engine`
    ///  (`Engine` is the storage scheme used by `GeMatrix`).  Hence, the
    ///  following would not compile:
    ///
    //  GeMatrixNoView M = FSView(4, 4, data, 4);
    //  error: This would not compile because
    //         GeMatrixNoView::Engine is of type "FullStorage<...>"
    //         FSView                 is of type "FullStorageView<...>"

    ///
    ///  Let us demonstrate what is a view and what not:
    ///    - We change `M(2,2)`,
    ///    - we change `B(1,2)`,
    ///    - we output `A`.
    ///  Only the change of `B(1,2)` affects elements of `A` as both matrices
    ///  share data.
    ///
    M(2,2) = -666;
    B(1,2) =  666;
    cout << "now: A = " << A << endl;
    cout << "Data dump of C-array: " << endl;
    for (int i=0; i<16; ++i) {
        cout << data[i] << "  ";
    }
    cout << endl;

    return 0;
}
