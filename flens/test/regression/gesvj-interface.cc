//
//  Regression test for lapack::svj via the Fortran interface dgesvj_ of
//  FLENS-LAPACK (flens/lapack/interface/src/gesvj.cc, part of
//  lapack_FLENS.a).
//
//  For triangular matrices with more than 64 columns svj uses a quasi block
//  cycling and created views of V for each block.  If V is not referenced
//  (JOBV='N', V has no rows) these views were out of range and svj aborted.
//  The reference DGEJSV calls DGESVJ like this, so the Netlib test suite
//  dsvd crashed as soon as DGEJSV was used with lapack_FLENS.a.
//
//  The singular values computed with JOBV='N' must agree with those for
//  JOBV='V' and with the singular values from DGESVD of the LAPACK library.
//
//  Requires a LAPACK library (build with REFLAPACK=...).  The dgesvj_
//  tested here is the one compiled from the FLENS sources included below.
//  As for lapack_FLENS.a, CBLAS is used (WITH_REFBLAS): FLENS-BLAS
//  implements rotm only through CBLAS.
//

#ifndef WITH_REFBLAS
#   define WITH_REFBLAS
#endif

#include <flens/lapack/interface/src/gesvj.cc>

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/iostream.h>

using namespace std;

extern "C" void
dgesvd_(const char *JOBU, const char *JOBVT, const int *M, const int *N,
        double *A, const int *LDA, double *S, double *U, const int *LDU,
        double *VT, const int *LDVT, double *WORK, const int *LWORK,
        int *INFO);

static int failures = 0;

static void
check(bool ok, const char *what, char typeA, char jobU, char jobV, double value)
{
    if (!ok) {
        ++failures;
        cerr << "FAILED: " << what << " (JOBA=" << typeA << ", JOBU=" << jobU
             << ", JOBV=" << jobV << "): " << value << endl;
    }
}

//
//  Singular values of a triangular n x n matrix computed by dgesvj_, sorted
//  in decreasing order.
//
static vector<double>
svj(char typeA, char jobU, char jobV, const vector<double> &A0, int n,
    int &info)
{
    vector<double> A = A0, sva(n), V(jobV=='V' ? n*n : 1);
    int ldA = n, mv = 0, ldV = (jobV=='V') ? n : 1;
    int lWork = max(6, 2*n);
    vector<double> work(lWork);

    flens::lapack::dgesvj_(&typeA, &jobU, &jobV, &n, &n, A.data(), &ldA,
                           sva.data(), &mv, V.data(), &ldV,
                           work.data(), &lWork, &info);
    for (int i=0; i<n; ++i) {
        sva[i] *= work[0];          // WORK(1) is the scaling factor
    }
    sort(sva.begin(), sva.end(), greater<double>());
    return sva;
}

static void
test(int n, char typeA, int seed)
{
    const double eps = numeric_limits<double>::epsilon();

    srand(seed);
    vector<double> A(n*n, 0.0);
    for (int j=0; j<n; ++j) {
        for (int i=0; i<n; ++i) {
            if ((typeA=='U' && i<=j) || (typeA=='L' && i>=j)) {
                A[i+j*n] = double(rand()%201-100)/10;
            }
        }
        A[j+j*n] += (j%2==0) ? 20 : -20;
    }

//  reference: singular values from DGESVD of the LAPACK library
    vector<double> B = A, s(n), work(1);
    int lWork = -1, info;
    double dummy;
    int one = 1;
    dgesvd_("N", "N", &n, &n, B.data(), &n, s.data(), &dummy, &one,
            &dummy, &one, &dummy, &lWork, &info);
    lWork = int(dummy);
    work.resize(lWork);
    dgesvd_("N", "N", &n, &n, B.data(), &n, s.data(), &dummy, &one,
            &dummy, &one, work.data(), &lWork, &info);
    check(info==0, "DGESVD: INFO", typeA, '-', '-', info);

    const char jobUs[2] = { 'U', 'N' };
    const char jobVs[2] = { 'N', 'V' };

    for (int u=0; u<2; ++u) {
        for (int v=0; v<2; ++v) {
            const char jobU = jobUs[u], jobV = jobVs[v];

            vector<double> sva = svj(typeA, jobU, jobV, A, n, info);
            check(info==0, "dgesvj_: INFO", typeA, jobU, jobV, info);

            double maxErr = 0;
            for (int i=0; i<n; ++i) {
                maxErr = max(maxErr, abs(sva[i]-s[i]));
            }
            const double tol = 10*n*eps*s[0];
            check(maxErr<=tol, "max |sigma - sigma(DGESVD)|",
                  typeA, jobU, jobV, maxErr);
        }
    }
}

int
main()
{
//  n > max(64, 4*kbl) selects the quasi block cycling for triangular input
    const int dims[3] = { 70, 100, 131 };

    for (int k=0; k<3; ++k) {
        test(dims[k], 'U', k+1);
        test(dims[k], 'L', k+11);
    }
//  small matrices use the plain one-sided Jacobi iteration
    test(10, 'U', 21);
    test(10, 'L', 22);

    if (failures) {
        cerr << failures << " failure(s)" << endl;
        return 1;
    }
    cout << "gesvj-interface: all tests passed" << endl;
    return 0;
}
