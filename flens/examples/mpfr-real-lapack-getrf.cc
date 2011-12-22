#include <iostream>

///
/// Include header file for `mpfr::real` before `flens.cxx` and make
/// sure that conversion operators are enabled
///
#define REAL_ENABLE_CONVERSION_OPERATORS
#include <external/real.hpp>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

///
///  Make typedef for using mpfr::real
///
typedef mpfr::real<53>   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const Underscore<IndexType> _;

    const IndexType m = 4,
                    n = 5;

    Matrix            Ab(m, n);
    IndexVector       piv(m);

    Ab =  2,   3,  -1,   0,  20,
         -6,  -5,   0,   2, -33,
          2,  -5,   6,  -6, -43,
          4,   6,   2,  -3,  49;

    cout << "Ab = " << Ab << endl;

    lapack::trf(Ab, piv);

    const auto A = Ab(_,_(1,m));
    auto       B = Ab(_,_(m+1,n));

    blas::sm(Left, NoTrans, T(1), A.upper(), B);

    cout << "X = " << B << endl;
}
