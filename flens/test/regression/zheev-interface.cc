//
//  Regression test for the Fortran interface zheev_ of FLENS-LAPACK
//  (flens/lapack/interface/src/heev.cc, part of lapack_FLENS.a).
//
//  For N==0 the interface created a view of length 3*N-2 == -2 on RWORK and
//  aborted (the Netlib test suite zsep crashed).  LAPACK documents the
//  length of RWORK as max(1, 3*N-2).
//
//  Requires a LAPACK library for the routines FLENS forwards to (build with
//  REFLAPACK=...).  The zheev_ tested here is the one compiled from the
//  FLENS sources included below.
//

#include <flens/lapack/interface/src/heev.cc>

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/iostream.h>

using namespace std;

typedef complex<double>  Z;

static int failures = 0;

static void
check(bool ok, const char *what, double value)
{
    if (!ok) {
        ++failures;
        cerr << "FAILED: " << what << ": " << value << endl;
    }
}

static void
testEmpty()
{
    const char jobz[2] = { 'N', 'V' };
    const char uplo[2] = { 'U', 'L' };

    for (int j=0; j<2; ++j) {
        for (int u=0; u<2; ++u) {
            int     n = 0, ldA = 1, lWork = 1, info = -1;
            Z       A[1], work[1];
            double  w[1], rWork[1];

            flens::lapack::zheev_(&jobz[j], &uplo[u], &n,
                                  reinterpret_cast<double *>(A), &ldA, w,
                                  reinterpret_cast<double *>(work), &lWork,
                                  rWork, &info);
            check(info==0, "zheev_ with N=0: INFO", info);
        }
    }
}

//
//  Small Hermitian matrices: eigenvalues of JOBZ='N' and JOBZ='V' agree and
//  the eigenvectors satisfy A*V = V*diag(w).
//
static void
testSmall(int n, int seed)
{
    const double eps = numeric_limits<double>::epsilon();

    srand(seed);
    vector<Z> A(n*n), V, B;
    for (int j=0; j<n; ++j) {
        for (int i=0; i<=j; ++i) {
            double re = double(rand()%201-100)/10;
            double im = (i==j) ? 0 : double(rand()%201-100)/10;
            A[i+j*n] = Z(re, im);
            A[j+i*n] = Z(re, -im);
        }
    }
    V = A;
    B = A;

    int ldA = n, lWork = 2*n, info;
    vector<Z>      work(lWork);
    vector<double> wV(n), wN(n), rWork(std::max(1, 3*n-2));

    flens::lapack::zheev_("V", "U", &n, reinterpret_cast<double *>(&V[0]),
                          &ldA, &wV[0], reinterpret_cast<double *>(&work[0]),
                          &lWork, &rWork[0], &info);
    check(info==0, "zheev_('V') INFO", info);
    flens::lapack::zheev_("N", "L", &n, reinterpret_cast<double *>(&B[0]),
                          &ldA, &wN[0], reinterpret_cast<double *>(&work[0]),
                          &lWork, &rWork[0], &info);
    check(info==0, "zheev_('N') INFO", info);

    double normA = 0;
    for (int k=0; k<n*n; ++k) {
        normA = std::max(normA, abs(A[k]));
    }
    normA *= n;
    const double bound = 64*n*eps*normA;

    for (int k=0; k<n; ++k) {
        check(abs(wV[k]-wN[k])<=bound, "eigenvalues of JOBZ='N' vs 'V'",
              abs(wV[k]-wN[k]));
    }
    double res = 0;
    for (int j=0; j<n; ++j) {
        for (int i=0; i<n; ++i) {
            Z av(0);
            for (int k=0; k<n; ++k) {
                av += A[i+k*n]*V[k+j*n];
            }
            res = std::max(res, abs(av - V[i+j*n]*wV[j]));
        }
    }
    check(res<=bound, "residual |A*V-V*W|", res);
}

int
main()
{
    testEmpty();
    for (int n=1; n<=12; ++n) {
        testSmall(n, 100+n);
    }

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "zheev-interface: all checks passed" << endl;
    return 0;
}
