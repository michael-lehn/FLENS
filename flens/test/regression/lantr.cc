//
//  Regression test for lapack::lan applied to triangular matrices.
//
//  For the infinity norm of an upper triangular matrix with non-unit
//  diagonal the inner loop incremented the column index j instead of the
//  row index i.  Only the first row was summed up, and the loop ran past
//  the last column of A.
//
//  All four norms are checked for upper/lower and unit/non-unit triangular
//  matrices against a straightforward computation
//  and, with -DWITH_REFERENCE_LAPACK, against dlantr/zlantr.
//

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
#include <cxxstd/vector.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

#ifdef WITH_REFERENCE_LAPACK
extern "C" {
double dlantr_(const char *norm, const char *uplo, const char *diag,
               const int *m, const int *n, const double *A, const int *ldA,
               double *work);
double zlantr_(const char *norm, const char *uplo, const char *diag,
               const int *m, const int *n, const complex<double> *A,
               const int *ldA, double *work);
}

double
lantr_(const char *norm, const char *uplo, const char *diag, const int *m,
       const int *n, const double *A, const int *ldA, double *work)
{
    return dlantr_(norm, uplo, diag, m, n, A, ldA, work);
}

double
lantr_(const char *norm, const char *uplo, const char *diag, const int *m,
       const int *n, const complex<double> *A, const int *ldA, double *work)
{
    return zlantr_(norm, uplo, diag, m, n, A, ldA, work);
}
#endif

static int failures = 0;

static void
check(bool ok, const char *what, double value, double expected)
{
    if (!ok) {
        ++failures;
        cerr << "FAILED: " << what << ": " << value
             << " (expected " << expected << ")" << endl;
    }
}

//
//  Straightforward computation of the norm of the triangular part of A.
//
template <typename MA>
double
naiveNorm(lapack::Norm norm, const MA &A, StorageUpLo upLo, Diag diag)
{
    const int m = A.numRows(), n = A.numCols();
    double maxAbs = 0, frob = 0, one = 0, inf = 0;
    vector<double> rowSum(m+1, 0.0), colSum(n+1, 0.0);

    for (int i=1; i<=m; ++i) {
        for (int j=1; j<=n; ++j) {
            bool inTriangle = (upLo==Upper) ? (i<=j) : (i>=j);
            if (!inTriangle) {
                continue;
            }
            double a = (i==j && diag==Unit) ? 1.0 : abs(A(i,j));
            maxAbs     = std::max(maxAbs, a);
            frob      += a*a;
            rowSum[i] += a;
            colSum[j] += a;
        }
    }
    for (int i=1; i<=m; ++i) {
        inf = std::max(inf, rowSum[i]);
    }
    for (int j=1; j<=n; ++j) {
        one = std::max(one, colSum[j]);
    }
    switch (norm) {
        case lapack::MaximumNorm:   return maxAbs;
        case lapack::OneNorm:       return one;
        case lapack::InfinityNorm:  return inf;
        default:                    return sqrt(frob);
    }
}

template <typename T>
static void
test(int m, int n)
{
    typedef GeMatrix<FullStorage<T> >   Matrix;
    typedef DenseVector<Array<double> > Vector;

    const lapack::Norm norms[4] = { lapack::MaximumNorm, lapack::OneNorm,
                                    lapack::InfinityNorm,
                                    lapack::FrobeniusNorm };
    const StorageUpLo  upLos[2] = { Upper, Lower };
    const Diag         diags[2] = { NonUnit, Unit };

    Matrix A(m, n);
    for (int i=1; i<=m; ++i) {
        for (int j=1; j<=n; ++j) {
            // row i == 2 dominates the infinity norm of the upper triangle
            A(i,j) = T(i==2 ? 50*j : i+j);
        }
    }

    for (int k=0; k<4; ++k) {
        for (int u=0; u<2; ++u) {
            for (int d=0; d<2; ++d) {
                Vector work(m);
                double result, expected;

                if (upLos[u]==Upper) {
                    auto U = (diags[d]==Unit) ? A.upperUnit() : A.upper();
                    result = lapack::lan(norms[k], U, work);
                } else {
                    auto L = (diags[d]==Unit) ? A.lowerUnit() : A.lower();
                    result = lapack::lan(norms[k], L, work);
                }
                expected = naiveNorm(norms[k], A, upLos[u], diags[d]);
                check(abs(result-expected)<=1e-13*expected,
                      "lan(TrMatrix) vs. naive computation",
                      result, expected);

#ifdef WITH_REFERENCE_LAPACK
                char normC = char(norms[k]);
                char upLoC = (upLos[u]==Upper) ? 'U' : 'L';
                char diagC = (diags[d]==Unit) ? 'U' : 'N';
                int  ldA   = A.leadingDimension();
                Vector rwork(m);
                double ref = lantr_(&normC, &upLoC, &diagC, &m, &n,
                                    A.data(), &ldA, rwork.data());
                check(abs(result-ref)<=1e-13*ref,
                      "lan(TrMatrix) vs. reference ?lantr", result, ref);
#endif
            }
        }
    }
}

int
main()
{
    for (int n=1; n<=7; ++n) {
        test<double>(n, n);
        test<complex<double> >(n, n);
    }

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "lantr: all checks passed" << endl;
    return 0;
}
