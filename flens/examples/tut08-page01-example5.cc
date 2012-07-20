#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef TpMatrix<PackedStorage<T, Upper> >  Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;

    ///
    /// The parameter define the dimensions (always a square matrix)
    ///

    Matrix         Ap(n);
    Vector         b(n), x(n);
    IndexVector    piv(n);
    
    for (IndexType i = 1; i <= n; ++i)
        for (IndexType j = i; j <= n; ++j)
            Ap(i,j) = i + 2*j;

    b = 20,
       -33,
       -43,
        49;

    cout << "Ap = " << Ap.symmetric() << endl;
    cout << "b  = " << b << endl;

    ///
    /// Compute $x = 0.5*(Ap+Ap^T) b$
    /// 
    x = Ap.symmetric()*b;
    
    cout << "x = " << x << endl;
    
    ///
    /// Solve $x = 0.5*(Ap+Ap^T) b$
    /// 
    lapack::sv(Ap.symmetric(), piv, x);

    cout << "b = " << b << endl;

    return 0;
}