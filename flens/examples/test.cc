#include <iostream>


/*
#define CXXBLAS_DEBUG_OUT(X)                                                \
   flens::verbose::ClosureLog::append(false) << "   => BLAS: " << X << ";"; \
*/


#include <flens/debug/aux/aux.h>
#include <flens/flens.cxx>
#include <flens/debug/aux/aux.tcc>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >   GeMatrix;
    typedef DenseVector<Array<T> >      DenseVector;
    typedef DenseVector::IndexType      IndexType;

    const Underscore<IndexType> _;

    int n = 3;

    verbose::ClosureLog::start("log");

    GeMatrix  A_(n,n), B_(n,n), C_(n,n), D_(n,n), E_(n,n);

    GeMatrix::View A = A_(_,_(1,3)),
                   B = B_(_,_(1,3)),
                   C = C_(_,_(1,3)),
                   D = D_(_,_(1,3)),
                   E = E_(_,_(1,3));

    A = 1;
    B = 2;
    C = 3;
    D = 4;
    E = 5;

    DenseVector    a(n), b(n), c(n), d(n), e(n);

    a = 1;
    b = 2;
    c = 3;
    d = 4;
    e = 5;

    /*
    b = A.upper()*b;
    c = A.upper()*b;
    c += A.upper()*b;

    C = transpose(A*B);
    C = A*B;
    A = A*A;
    C = C + transpose(A)*B;
    C = 2*C + transpose(A)*B;
    C = 2*C + transpose(transpose(A)*B);
    C = transpose(2*C + transpose(A)*B);

    b = 2*b + 3*transpose(A)*a;
    b = 3*transpose(A)*a + 2*b;
    */

    cout << "HasEngine<ScalarValue<int> >::value = "
         << HasEngine<ScalarValue<int> >::value
         << endl;

    cout << "HasEngine<GeMatrix>::value = "
         << HasEngine<GeMatrix>::value
         << endl;

    cout << "HasEngine<GeMatrix::TrMatrixView>::value = "
         << HasEngine<GeMatrix::TriangularView>::value
         << endl;

    cout << "IsFullStorage<GeMatrix::Engine>::value = "
         << IsFullStorage<GeMatrix::Engine>::value
         << endl;

    cout << "HasFullStorage<GeMatrix>::value = "
         << HasFullStorage<GeMatrix>::value
         << endl;

    cout << "IsFullStorage<int>::value = "
         << IsFullStorage<int>::value
         << endl;

    using DEBUGCLOSURE::identical;
    cout << "identical(C, C.upper()) = " << identical(C, C.upper()) << endl;

    C = A*C.upper();
    // C = C.upper()*A;
    // C = A.upper()*C.upper();

    /*
    C = A.upper();

    a = c - d;

    c = b - A*a;
    c = b - 2*A*a;
    c = 2*c + transpose(2*transpose(A))*a;
    c = 2*c + transpose(2*transpose(A+B))*(a+b);

    e = 2*b + A*(a+c);
    cout << "e = " << e << endl;

    e = 2*b + transpose(A)*(a+c);
    cout << "e = " << e << endl;

    e = 2*b + transpose(3*A)*(a+c);
    cout << "e = " << e << endl;

    e = 2*b + transpose(3*A + transpose(A))*(a+c);
    cout << "e = " << e << endl;


    cout << "A = " << A << endl;
    A += transpose(A);
    cout << "A = " << A << endl;

    A += transpose(A);
    cout << "A = " << A << endl;

    GeMatrix X(n, n);
    X = A;
    cout << "X = " << X << endl;

    X(_(1,2),_).upper() = transpose(E(_,_(1,2)).lower());
    cout << "X = " << X << endl;


    b = A*b;

    b = A.upper()*b;
    c += 3*(2*transpose(A.upper())*b);
    b += A.upper()*b;

    B = C.lower();
    cout << "A = " << A << endl;
    cout << "C = " << C << endl;
    cout << "A.upper() = " << A.upper() << endl;
    cout << "C.lower() = " << C.lower() << endl;
    cout << "B = " << B << endl;

    D = C.lower();
    cout << "D = "  << D << endl;


    A = 3*B + (4*C + 5*D);
    cerr << "A = " << A << endl;

    A = 3*B + (4*C + 2*(5*D));
    cerr << "A = " << A << endl;

    A = 3*B + (4*C + D/5);
    cerr << "A = " << A << endl;

    A = 3*B + (4*C + D/5/2);
    cerr << "A = " << A << endl;

    A = B + (C + D);
    cerr << "A = " << A << endl;

    A = B + (C - D);
    cerr << "A = " << A << endl;

    B(1,_) = 4;
    cerr << "B = " << B << endl;

    A = transpose(B);
    cerr << "A = " << A << endl;

    A = transpose(A);
    cerr << "A = " << A << endl;

    A = transpose(A(_(1,n),_(1,n)));
    cerr << "A = " << A << endl;

    A = 3*B;
    cerr << "A = " << A << endl;

    A = 3*A;
    cerr << "A = " << A << endl;

    A = 2*(3*B);
    cerr << "A = " << A << endl;

    A = B/3;
    cerr << "A = " << A << endl;

    A = A/3;
    cerr << "A = " << A << endl;

    A = B/3/2;
    cerr << "A = " << A << endl;


    a = a + b;

    a = b + a;

    a = a + c + d + e+ b + a(_(1,n));

    a = a + a(_(1,n)) + d + e+ b + a(_(1,n));

    a = 3*a;

    a = a*3*2;

    e = 2*a - b/2 + 2*c + 5*d;
    cerr << "e = " << e << endl;

    e = M_PI*(2*a) - b/2 + 2*(3*c) + d;
    cerr << "e = " << e << endl;

    e = M_PI*(2*a) - b/2/3 + 2*(3*c) + d;
    cerr << "e = " << e << endl;

    e = b/2;
    cerr << "e = " << e << endl;

    const DenseVector &x = a;
    cerr << "x = " << x << endl;

    const DenseVector &y = b;
    cerr << "y = " << y << endl;

    e = 2*x*y*x;
    cerr << "e = " << e << endl;

    e += (x+x)*(y*x);
    cerr << "e = " << e << endl;

    e = ((a - b) - c) - d;
    cerr << "e = " << e << endl;

    e = a - (b - (c + d));
    cerr << "e = " << e << endl;

    e = a + (b + (c + d));
    cerr << "e = " << e << endl;


    D = transpose(C.lower())*transpose(A+B);

    cout << "D = "  << D << endl;
    cout << "C = "  << C << endl;
    cout << "C.lower().symmetric() = " << C.lower().symmetric() << endl;
    GeMatrix XX = C.lower().symmetric() + C;
    cout << "-> XX = "  << XX << endl;

    D = C*A;
    cout << "D = "  << D << endl;
    D = C.lower().symmetric()*D;
    cout << "D = "  << D << endl;
    D = D.lower().symmetric()*A;

    cout << "D = "  << D << endl;
    b = (D.lower().symmetric()+D)*a;
    b = D.lower().symmetric()*a;

    A = E.lower() + C.upper();
    cout << "E.lower() = "  << E.lower() << endl;
    cout << "C.upper() = "  << C.upper() << endl;
    cout << "A = "  << A << endl;
    */

    verbose::ClosureLog::stop();
}
