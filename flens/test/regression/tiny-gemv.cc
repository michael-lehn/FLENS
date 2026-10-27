//
//  Regression test for cxxblas::gemv on tiny matrices.
//
//  1) For beta==0 the kernels computed y *= beta before accumulating A*x.
//     As in the reference BLAS, y must not be read in this case: it may be
//     uninitialized (TinyArray does not initialize its elements) or contain
//     NaN/Inf.  examples/tiny_gemv printed "nan" when compiled with clang.
//  2) The transposed variants advanced the index of x by incY instead of
//     incX.
//
//  Reference values are computed with straightforward loops.
//

#include <cxxstd/cmath.h>
#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
#include <cxxstd/limits.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

static int failures = 0;

static void
check(bool ok, const char *what)
{
    if (!ok) {
        ++failures;
        cerr << "FAILED: " << what << endl;
    }
}

//
//  The example from flens/examples/tiny_gemv.cc, with y set to NaN first.
//
static void
testTinyMatrixVector()
{
    typedef GeTinyMatrix<TinyFullStorage<double, 4, 4> >  Matrix;
    typedef TinyVector<TinyArray<double, 4> >             Vector;

    Matrix A;
    A(1,1) = 1;  A(1,2) = 2;  A(1,3) = 3;  A(1,4) =  4;
    A(2,1) = 5;  A(2,2) = 6;  A(2,3) = 7;  A(2,4) =  8;
    A(3,1) = 9;  A(3,2) = 8;  A(3,3) = 7;  A(3,4) =  6;
    A(4,1) = 5;  A(4,2) = 4;  A(4,3) = 3;  A(4,4) = 20;

    Vector v, w;
    v.fill(1);

    w.fill(numeric_limits<double>::quiet_NaN());
    w = A*v;
    const double Av[4] = { 10, 26, 30, 32 };
    for (int i=1; i<=4; ++i) {
        check(w(i)==Av[i-1], "tiny matrix: w = A*v");
    }

    w.fill(numeric_limits<double>::quiet_NaN());
    w = transpose(A)*v;
    const double ATv[4] = { 20, 20, 20, 38 };
    for (int i=1; i<=4; ++i) {
        check(w(i)==ATv[i-1], "tiny matrix: w = A^T*v");
    }
}

//
//  The kernel itself with incX = 2, incY = 3 for all transpose variants and
//  beta = 0 (y initialized with NaN) as well as beta != 0.
//
template <typename T>
static void
testKernel(Transpose trans, T beta)
{
    const int m = 2, n = 3, ldA = 4;
    const int incX = 2, incY = 3;
    const bool  transA = (trans==Trans || trans==ConjTrans);
    const bool  conjA  = (trans==Conj || trans==ConjTrans);
    const int   lenX   = transA ? m : n;
    const int   lenY   = transA ? n : m;

    // row major m x n matrix with leading dimension ldA
    T A[m*ldA];
    for (int i=0; i<m; ++i) {
        for (int j=0; j<ldA; ++j) {
            A[i*ldA+j] = (j<n) ? T(i+1, j-i) : T(-1000);
        }
    }

    T x[incX*3], y[incY*3], y0[incY*3];
    for (int k=0; k<incX*3; ++k) {
        x[k] = T(k+1, 2*k-1);
    }
    for (int k=0; k<incY*3; ++k) {
        y0[k] = (beta==T(0)) ? T(numeric_limits<double>::quiet_NaN())
                             : T(k-2, k);
        y[k]  = y0[k];
    }

    const T alpha(2, 1);
    cxxblas::gemv<m, n, T, ldA, T, incX, T, incY>(trans, alpha, A, x,
                                                  beta, y);

    for (int k=0; k<lenY; ++k) {
        T sum = T(0);
        for (int l=0; l<lenX; ++l) {
            T a = transA ? A[l*ldA+k] : A[k*ldA+l];
            sum += (conjA ? conj(a) : a) * x[l*incX];
        }
        T expected = alpha*sum;
        if (beta!=T(0)) {
            expected += beta*y0[k*incY];
        }
        check(abs(y[k*incY]-expected)<=1e-13*abs(expected),
              "tiny gemv kernel: y = beta*y + alpha*op(A)*x");
    }
    // elements of y between the strided entries must not be touched
    for (int k=0; k<incY*3; ++k) {
        if (k%incY==0 && k/incY<lenY) {
            continue;
        }
        bool same = (y[k]==y0[k]) || (y[k]!=y[k] && y0[k]!=y0[k]);
        check(same, "tiny gemv kernel: no write outside of y");
    }
}

int
main()
{
    typedef complex<double> Z;

    testTinyMatrixVector();

    const Transpose trans[4] = { NoTrans, Trans, Conj, ConjTrans };
    for (int t=0; t<4; ++t) {
        testKernel<Z>(trans[t], Z(0));
        testKernel<Z>(trans[t], Z(0.5, -1));
    }

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "tiny-gemv: all checks passed" << endl;
    return 0;
}
