#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    GeMatrix<FullStorage<double> >   AB(4, 4+3);
    DenseVector<Array<int> >         piv(4);

    Underscore<int> _;

    auto A = AB(_,_(1,4));
    auto B = AB(_,_(5,7));

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

    fillRandom(B);

    cout << "(A,B) = " << AB << endl;

///
/// *A pripori* compute the triangular factorization $PLU = (A,B)$ using the
/// __FLENS-LAPACK__ function `lapack::trf`.
///
    lapack::trf(AB, piv);

///
/// Solve $U X = B$ using __FLENS-BLAS__ function `blas::sv.  Matrix $B$ gets
/// overwritten with the solution $X$.
///
    blas::sm(Left, NoTrans, 1.0, A.upper(), B);

    cout << "X = " << B << endl;
}

///
///     :links: FLENS-LAPACK -> doc:flens/lapack/lapack
///             FLENS-BLAS   -> doc:flens/blas/blas
///
