//
//  Regression test for cxxblas::copy and cxxblas::ccopy on tiny vectors.
//
//  Both used y[iX] instead of y[iY], i.e. they wrote with the stride of x.
//  For incX!=incY this gave wrong results and wrote out of bounds.  The
//  transposed copy of tiny matrices (B = transpose(A)) is implemented in
//  terms of copy with incX = leading dimension of A and incY = 1, so it was
//  affected as well.
//

#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
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

template <typename T>
static void
testKernel(bool conjugate)
{
    const int n = 3;
    // x has stride 2: {1, 2, 3} interleaved with -1
    T x[2*n] = { T(1,1), T(-1), T(2,2), T(-1), T(3,3), T(-1) };

    // y has stride 1, guarded by sentinels to detect out-of-bounds writes
    const T sentinel(-99);
    T buffer[n+2*n];
    for (int i=0; i<3*n; ++i) {
        buffer[i] = sentinel;
    }
    T *y = buffer + 1;

    if (conjugate) {
        cxxblas::ccopy<n, T, 2, T, 1>(x, y);
    } else {
        cxxblas::copy<n, T, 2, T, 1>(x, y);
    }

    for (int i=0; i<n; ++i) {
        T expected = conjugate ? conj(x[2*i]) : x[2*i];
        check(y[i]==expected, conjugate ? "ccopy: y(i) == conj(x(i))"
                                        : "copy: y(i) == x(i)");
    }
    check(buffer[0]==sentinel, "copy: no write before y");
    for (int i=n+1; i<3*n; ++i) {
        check(buffer[i]==sentinel, "copy: no write after y");
    }
}

static void
testTransposedTinyMatrix()
{
    typedef GeTinyMatrix<TinyFullStorage<double, 3, 4> >  Matrix34;
    typedef GeTinyMatrix<TinyFullStorage<double, 4, 3> >  Matrix43;

    Matrix34 A;
    for (int i=1; i<=3; ++i) {
        for (int j=1; j<=4; ++j) {
            A(i,j) = 10*i + j;
        }
    }
    Matrix43 B;
    B = transpose(A);

    for (int i=1; i<=3; ++i) {
        for (int j=1; j<=4; ++j) {
            check(B(j,i)==A(i,j), "tiny matrix: B = transpose(A)");
        }
    }
}

int
main()
{
    testKernel<complex<double> >(false);
    testKernel<complex<double> >(true);
    testTransposedTinyMatrix();

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "tiny-copy: all checks passed" << endl;
    return 0;
}
