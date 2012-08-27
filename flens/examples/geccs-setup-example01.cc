#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
///
/// The following typedef `Coord` is a shortcut for a coordiante storage formate
/// with the following properties:
///  - Element are of type `double`
///  - Indices
///     - are of type `int` and
///     - start at zero (like in C arrays).
///  - When elements get accumulated they get sorted *column-wise*.
///
    typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<double, CoordColRowCmp, IndexBase>  Coord;

///
/// We define a general sparse $m \times n$ matrix $A$.
///
    const IndexType m = 5;
    const IndexType n = 5;
    GeCoordMatrix<Coord>  A(m, n);

///
/// We add some values to some entries $a_{i_k, j_k}$ of $A$.  Arrows indicate
/// that certain entries occur more then once.
///
    A(0,0) += 2;
    A(0,1) += 3;
    A(1,0) += 3;
    A(1,2) += 4;
    A(1,4) += 3;       // <-
    A(1,4) += 3;       // <-
    A(2,1) -= 1;
    A(2,2) -= 3;
    A(2,3) += 2;
    A(3,2) += 1;
    A(4,1) += 2;       // <-
    A(4,2) += 2;
    A(4,4) += 1;
    A(4,1) += 2;       // <-

///
/// Next we convert the sparse matrix $A$ with *coordinate storage* to a sparse
/// matrix $B$ with *compressed column storage*.  Values in the storage format
/// of $A$ get accumulated
///
    GeCCSMatrix<CCS<double, IndexBase> >  B = A;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

///
/// You can continue to add values to entries in $A$.
///
    A(0,2) += 2;
    A(0,4) += 3;
    A(1,0) += 3;
    A(1,2) += 4;
    A(1,4) += 6;
    A(2,1) -= 1;
    A(2,2) -= 3;
    A(2,3) += 2;
    A(3,2) += 1;
    A(4,1) += 4;
    A(4,2) += 2;
    A(4,4) += 1;

///
/// Again we convert $A$ to *compressed column storage*.
///
    B = A;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

///
/// For debugging it's sometime convenient to densify the sparse matrix.  I.e.
/// we convert it to GeMatrix.
///
    GeMatrix<FullStorage<double> >  C = B;
    cout << "C = " << C << endl;
}
