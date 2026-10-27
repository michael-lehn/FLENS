//
//  Regression test for lapack::lqf with less than optimal workspace.
//
//  In this case lqf queries ilaenv for the minimal block size and passed a
//  null pointer instead of "" as the 'opts' string.  If FLENS is configured
//  to use the native ilaenv (as lapack_FLENS.a is: USE_ILAENV_WITH_UNDERSCORE,
//  USE_NATIVE_ILAENV), the string is passed on to the Fortran ilaenv and
//  strlen(nullptr) crashed.  The Netlib test suite zsvd crashed in zgelqf_
//  because of this.
//
//  The factorization computed with insufficient workspace is compared with
//  the one computed with optimal workspace and with the reference
//  LAPACK routines dgelqf/zgelqf.  Requires a LAPACK library (REFLAPACK).
//

#define USE_CXXLAPACK
#define USE_ILAENV_WITH_UNDERSCORE
#define USE_NATIVE_ILAENV

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

typedef complex<double>  Z;

static int failures = 0;

static void
check(bool ok, const char *what, double value, double bound)
{
    if (!ok) {
        ++failures;
        cerr << "FAILED: " << what << ": " << value
             << " (bound " << bound << ")" << endl;
    }
}

template <typename MA, typename MB>
static double
maxAbsDiff(const MA &A, const MB &B)
{
    double d = 0;
    for (int j=1; j<=A.numCols(); ++j) {
        for (int i=1; i<=A.numRows(); ++i) {
            d = std::max(d, double(abs(A(i,j)-B(i,j))));
        }
    }
    return d;
}

template <typename T>
static void
test(int m, int n, int lWork)
{
    typedef GeMatrix<FullStorage<T> >  Matrix;
    typedef DenseVector<Array<T> >     Vector;

    const double eps = numeric_limits<double>::epsilon();
    const int    k   = std::min(m, n);

    srand(m+n);
    Matrix A(m, n);
    for (int j=1; j<=n; ++j) {
        for (int i=1; i<=m; ++i) {
            A(i,j) = T(double(rand()%2001-1000)/1000);
        }
    }
    Matrix A1 = A, A2 = A, A3 = A;
    Vector tau1(k), tau2(k), tau3(k);

    // insufficient workspace: triggers the ilaenv query for NBMIN
    Vector work1(lWork);
    lapack::lqf(A1, tau1, work1);

    // optimal workspace (computed by lqf itself)
    Vector work2;
    lapack::lqf(A2, tau2, work2);

    // reference LAPACK with the same workspace as in the first call
    Vector work3(lWork);
    int info = cxxlapack::gelqf<int>(m, n, A3.data(), m, tau3.data(),
                                     work3.data(), lWork);
    check(info==0, "reference ?gelqf INFO", info, 0);

    const double bound = 64*std::max(m, n)*eps*std::max(m, n);
    check(maxAbsDiff(A1, A3)<=bound,
          "lqf (small workspace) vs. reference ?gelqf", maxAbsDiff(A1, A3),
          bound);
    check(maxAbsDiff(A1, A2)<=bound,
          "lqf (small workspace) vs. lqf (optimal workspace)",
          maxAbsDiff(A1, A2), bound);
    double dTau = 0;
    for (int i=1; i<=k; ++i) {
        dTau = std::max(dTau, double(abs(tau1(i)-tau3(i))));
    }
    check(dTau<=bound, "tau: lqf vs. reference ?gelqf", dTau, bound);
}

int
main()
{
    // k = min(m,n) must exceed the crossover point NX (128 for GELQF) and
    // lWork must be smaller than m*NB
    test<double>(150, 160, 2*150);
    test<double>(200, 140, 3*200);
    test<Z>(150, 160, 2*150);
    test<Z>(200, 140, 3*200);

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "lqf-ilaenv: all checks passed" << endl;
    return 0;
}
