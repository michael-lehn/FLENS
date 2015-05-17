#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    GeMatrix<FullStorage<double> >   A(4, 4);
    DenseVector<Array<double> >      b(4);
    DenseVector<Array<int> >         piv(4);

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

    cout << "A = " << A << endl;

///
/// *A pripori* compute the triangular factorization $PLU = A$ using the
/// __FLENS-LAPACK__ function `lapack::trf`.
///
    lapack::trf(A, piv);

///
/// The right hand sides are only known a posteriori.  Furthermore we have
/// three different right hand sides for the same coefficient matrix $A$.
///
/// *A posteriori* call the triangular solvers to solve $Ax=b$ for different
/// right hand sides $b$.
///
    for (int i=0; i<3; ++i) {
        cout << "Right hand side b in problem set " << i << endl;
        fillRandom(b);
        cout << "b = " << b << endl;
///
///     Apply the permutation: $b \leftarrow Pb$ using the __FLENS-LAPACK__
///     function `lapack::laswp`.
///
        lapack::laswp(b, piv);

///
///     Solve $L y = b$ using __FLENS-BLAS__ function `blas::sv`.  Vector $b$
///     will be overwritten with $y$.
///
        blas::sv(NoTrans, A.lowerUnit(), b);

///
///     Solve $U x = b$ using __FLENS-BLAS__ function `blas::sv`.  Vector $b$
///     will be overwritten with $x$.
///
        blas::sv(NoTrans, A.upper(), b);

        cout << "x = " << b << endl;
    }
}

///
///     :links: FLENS-LAPACK -> doc:flens/lapack/lapack
///             FLENS-BLAS   -> doc:flens/blas/blas
///
