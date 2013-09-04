#define USE_PLAYGROUND
#define WITH_MKLBLAS
#include <flens/flens.cxx>
#include <iostream>

using namespace std;
using namespace flens;

typedef complex<float> T;
int main(void)
{
	///
	/// The following typedef `Coord` is a shortcut for a coordiante storage formate
	/// with the following properties:
	///  - Element are of type `complex<float>`
	///  - Indices
	///     - are of type `int` and
	///     - start at one.
	///  - When elements get accumulated they get sorted *column-wise*.
	///
    typedef int                              IndexType;
    typedef CoordStorage<T, CoordRowColCmp>  Coord;


    ///
    /// We define a general sparse $n \times n$ matrix $A$
    /// and two dense vectors $x$ and $b$
    ///
	IndexType n = 5;
    DenseVector<Array<T> >   x(n), b(n);
    GeCoordMatrix<Coord>     A(n, n);


    ///
    /// We add some values to some entries $a_{i_k, j_k}$ of $A$.  Arrows indicate
    /// that certain entries occur more then once.
    ///

    A(1, 1) += T(3, 1);
    A(2, 1) += T(4, 2);
    A(4, 1) += T(6, -3);
    A(2, 2) += T(-3, 4);
    A(3, 2) += T(2, 1);
    A(3, 3) += T(1, 0);
    A(3, 4) += T(4, 2);
    A(1, 4) += T(-2, 3);
    A(3, 4) += T(1, 2);
    A(4, 4) += T(1, -1);
    A(5, 5) += T(2, 4);

    ///
    /// Next we convert the sparse matrix $A$ with *coordinate storage* to a sparse
    /// matrix $B$ with *compressed row storage*.
    ///
    GeCRSMatrix<CRS<T> >               B = A;

	///
    /// Setup the right-hand side
    ///
    b = T(3, 1), T(4, 1), T(5, 9), T(2, 6), T(5, 3);

	cout << b << endl;
    ///
    /// Solve $B \cdot x = b$
    ///
    mkldss::sv(ConjTrans, B, x, b);

    ///
    /// Check solution
    ///
	b = conjTrans(B)*x;
	cout << b << endl;

	return 0;
}
