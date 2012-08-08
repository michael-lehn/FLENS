#include <iostream>
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

int
main()
{
    typedef complex<double>                  Complex;
    const Complex                            I(0,1);

    typedef GeMatrix<FullStorage<Complex> >  GeMatrix;
    typedef typename GeMatrix::IndexType     IndexType;
    typedef DenseVector<Array<IndexType> >   IndexVector;
    typedef DenseVector<Array<Complex> >     DenseVector;

    GeMatrix        A(4,4);
    IndexVector     jPiv;
    DenseVector     tau;
    DenseVector     b(4);
    //DenseVector   work;

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

    b = 20,
       -33,
       -43,
        49;

///
/// Make that $A$ and $b$ somehow look complex ;-)
///
    A *= I;
    b *= 2.0*I;

    cout << "A = " << A << endl;
    cout << "b = " << b << endl;

///
/// Compute the factorization $AP = QR$.  As vector `jPiv` has initially length
/// zero it gets internally resized to length $n$ and initialized with Zero.
/// This means that *all columns are free columns*.  Also note that the
/// workspace gets created implicitly and temporarily.  So you might not want to do this inside a loop.
///
    lapack::qp3(A, jPiv, tau);
    //lapack::qp3(A, jPiv, tau, work);

///
/// Compute $\tilde{b} = Q^H b$. Vector $b$ gets overwritten with $\tilde{b}$.
///
    lapack::unmqr(Left, ConjTrans, A, tau, b);
    //lapack::unmqr(Left, ConjTrans, A, tau, b, work);

///
/// Solve $R u = \tilde{b}$.  Vector $b$ gets overwritten with $u =P^T x$.
///
    blas::sv(NoTrans, A.upper(), b);

///
/// Compute $x = Pu$.  Note that we need an extra vector here!
///
    DenseVector     x(4);
    for (IndexType i=1; i<=x.length(); ++i) {
        x(jPiv(i)) = b(i);
    }

    cout << "x = " << x << endl;

///
/// Output some additional information:
///
    cout << "jPiv = " << jPiv << endl;
}
