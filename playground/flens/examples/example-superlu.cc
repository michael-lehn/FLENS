#include <cxxstd/iostream.h>

#define USE_PLAYGROUND
#define WITH_SUPERLU
#include <flens/flens.cxx>

using namespace std;
using namespace flens;


typedef double T;

int
main()
{
    ///
    /// The following typedef `Coord` is a shortcut for a coordiante storage
    /// formate with the following properties:
    ///  - Element are of type `double`
    ///  - Indices
    ///     - are of type `int` and
    ///     - start at zero (like in C arrays).
    ///  - When elements get accumulated they get sorted *column-wise*.
    ///
    typedef int                                              IndexType;
    typedef IndexBaseZero<IndexType>                         IndexBase;
    typedef CoordStorage<T, CoordColRowCmp, IndexBase>  Coord;


    ///
    /// We define a general sparse $n \times n$ matrix $A$
    /// and two dense vectors $x$ and $b$
    ///
    IndexType n = 5;
    DenseVector<Array<T, IndexBase> >         x(n), b(n);
    DenseVector<Array<IndexType, IndexBase> > pc, pr;
    GeCoordMatrix<Coord>              A(n, n);


    ///
    /// We add some values to some entries $a_{i_k, j_k}$ of $A$.  Arrows
    /// indicate that certain entries occur more then once.
    ///

    A(0, 0) += 2;
    A(0, 1) += 3;
    A(1, 0) += 3;
    A(1, 2) += 4;
    A(1, 4) += 6;
    A(2, 1) += -1;
    A(2, 2) += -3;
    A(2, 3) += 2;
    A(3, 2) += 1;
    A(4, 1) += 4;
    A(4, 2) += 2;
    A(4, 4) += 1;

    ///
    /// Next we convert the sparse matrix $A$ with *coordinate storage* to a
    /// sparse matrix $B$ with *compressed column storage*.
    ///
    GeCCSMatrix<CCS<T> >               B = A;

    ///
    /// Setup the right-hand side
    ///
    x = 1, 2, 3, 4, 5;
    cout << x << endl;

    ///
    /// Solve $B \cdot x = b$
    ///
    superlu::sv(B, pc, pr, x);

    ///
    /// Check the solution
    ///

    b = B*x;
    cout << b << endl;

    return 0;
}
