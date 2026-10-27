//
//  Regression test for lapack::sterf and the eigenvalues-only path of
//  lapack::ev for Hermitian matrices (which calls sterf).
//
//  sterf used to return the unreduced diagonal: the loop searching for a
//  small subdiagonal element in the QL branch declared a new loop variable
//  and incremented the wrong one.
//
//  Also, the bound for scaling was computed as One/SafeMin (the enumerator
//  of MachineParameter) instead of One/safeMin.  Blocks with large norm
//  were therefore scaled down so much that the remaining entries became
//  subnormal, see testWideRange().
//
//  Reference values:
//    - Symmetric tridiagonal Toeplitz matrices have the eigenvalues
//      a + 2*b*cos(k*pi/(n+1)), k=1..n.
//    - For Hermitian matrices, ev(false, ...) must agree with ev(true, ...)
//      and the eigenvectors computed by ev(true, ...) must satisfy
//      A*V = V*diag(w) and V^H*V = I.
//    - With -DWITH_REFERENCE_LAPACK the results are also compared with the
//      external LAPACK routines dsterf and zheev.
//

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/complex.h>
#include <cxxstd/cstdlib.h>
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

typedef double                              T;
typedef complex<T>                          Z;
typedef DenseVector<Array<T> >              DVector;
typedef DenseVector<Array<Z> >              ZVector;
typedef GeMatrix<FullStorage<Z> >           ZGeMatrix;
typedef DVector::IndexType                  Index;

#ifdef WITH_REFERENCE_LAPACK
extern "C" {
void dsterf_(const int *n, double *d, double *e, int *info);
void zheev_(const char *jobz, const char *uplo, const int *n, Z *A,
            const int *ldA, double *w, Z *work, const int *lWork,
            double *rWork, int *info);
}
#endif

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

static T
maxAbsDiff(const DVector &x, const DVector &y)
{
    T d = 0;
    for (Index i=1; i<=x.length(); ++i) {
        d = std::max(d, std::abs(x(i)-y(i)));
    }
    return d;
}

//
//  sterf on tridiagonal Toeplitz matrices with known spectrum
//
static void
testToeplitz(Index n, T a, T b)
{
    const T eps = numeric_limits<T>::epsilon();
    const T pi  = std::acos(T(-1));

    DVector d(n), e(n-1), exact(n);
    for (Index i=1; i<=n; ++i) {
        d(i) = a;
    }
    for (Index i=1; i<=n-1; ++i) {
        e(i) = b;
    }
    for (Index k=1; k<=n; ++k) {
        exact(k) = a + 2*b*std::cos(k*pi/(n+1));
    }
    std::sort(exact.data(), exact.data()+n);    // sterf sorts ascending

    DVector d0 = d, e0 = e;
    Index info = lapack::sterf(d, e);

    const T bound = 32*n*eps*(std::abs(a)+2*std::abs(b));
    check(info==0, "sterf info", info, 0);
    check(maxAbsDiff(d, exact)<=bound, "sterf vs. exact eigenvalues",
          maxAbsDiff(d, exact), bound);

#ifdef WITH_REFERENCE_LAPACK
    int n_ = n, info_;
    dsterf_(&n_, d0.data(), e0.data(), &info_);
    check(info_==0, "dsterf info", info_, 0);
    check(maxAbsDiff(d, d0)<=bound, "sterf vs. reference dsterf",
          maxAbsDiff(d, d0), bound);
#endif
}

//
//  One huge diagonal entry D weakly coupled to the 2x2 block [1 1; 1 1].
//  The coupling e(1) = sqrt(D)*1e-10 is large enough to prevent a split but
//  only perturbs the eigenvalues by about e(1)^2/D = 1e-20.  Hence the
//  eigenvalues are 0, 2 and D (up to rounding).
//
static void
testWideRange(T D)
{
    const T eps = numeric_limits<T>::epsilon();

    DVector d(3), e(2);
    d = D, 1, 1;
    e = std::sqrt(D)*1e-10, 1;

#ifdef WITH_REFERENCE_LAPACK
    DVector d0 = d, e0 = e;
#endif

    Index info = lapack::sterf(d, e);
    check(info==0, "sterf info (wide range)", info, 0);
    check(std::abs(d(1))<=16*eps, "sterf wide range: eigenvalue 0",
          std::abs(d(1)), 16*eps);
    check(std::abs(d(2)-2)<=16*eps, "sterf wide range: eigenvalue 2",
          std::abs(d(2)-2), 16*eps);
    check(std::abs(d(3)-D)<=16*eps*D, "sterf wide range: eigenvalue D",
          std::abs(d(3)-D), 16*eps*D);

#ifdef WITH_REFERENCE_LAPACK
    int n_ = 3, info_;
    dsterf_(&n_, d0.data(), e0.data(), &info_);
    check(info_==0, "dsterf info (wide range)", info_, 0);
    for (Index i=1; i<=2; ++i) {
        check(std::abs(d(i)-d0(i))<=16*eps,
              "sterf vs. reference dsterf (wide range)",
              std::abs(d(i)-d0(i)), 16*eps);
    }
    check(std::abs(d(3)-d0(3))<=16*eps*D,
          "sterf vs. reference dsterf (wide range, D)",
          std::abs(d(3)-d0(3)), 16*eps*D);
#endif
}

//
//  Hermitian eigenvalue problem: eigenvalues-only vs. eigenvectors mode
//
static void
testHermitian(Index n, int seed)
{
    const T eps = numeric_limits<T>::epsilon();

    srand(seed);
    ZGeMatrix A(n, n);
    for (Index j=1; j<=n; ++j) {
        for (Index i=1; i<=j; ++i) {
            T re = T(rand()%2001 - 1000)/100;
            T im = (i==j) ? T(0) : T(rand()%2001 - 1000)/100;
            A(i,j) = Z(re, im);
            A(j,i) = Z(re, -im);
        }
    }
    T normA = 0;
    for (Index j=1; j<=n; ++j) {
        for (Index i=1; i<=n; ++i) {
            normA = std::max(normA, std::abs(A(i,j)));
        }
    }
    normA *= n;

    ZGeMatrix V = A, B = A;
    DVector      wV(n), wN(n);
    ZVector   work;
    DVector      rWork;

    Index infoV = lapack::ev(true, V.upper().hermitian(), wV, work, rWork);
    Index infoN = lapack::ev(false, B.upper().hermitian(), wN, work, rWork);
    check(infoV==0, "ev(true) info", infoV, 0);
    check(infoN==0, "ev(false) info", infoN, 0);

    const T bound = 64*n*eps*normA;
    check(maxAbsDiff(wV, wN)<=bound, "ev(false) vs. ev(true) eigenvalues",
          maxAbsDiff(wV, wN), bound);

    // eigenvalues must be in ascending order
    for (Index i=1; i<n; ++i) {
        check(wN(i)<=wN(i+1), "ev(false) ascending order", wN(i), wN(i+1));
    }

    // residual ||A*V - V*diag(w)|| and orthogonality ||V^H*V - I||
    T res = 0, orth = 0;
    for (Index j=1; j<=n; ++j) {
        for (Index i=1; i<=n; ++i) {
            Z av(0), vhv(0);
            for (Index k=1; k<=n; ++k) {
                av  += A(i,k)*V(k,j);
                vhv += conj(V(k,i))*V(k,j);
            }
            res  = std::max(res, std::abs(av - V(i,j)*wV(j)));
            orth = std::max(orth, std::abs(vhv - Z(i==j ? 1 : 0)));
        }
    }
    check(res<=bound, "ev(true) residual |A*V-V*W|", res, bound);
    check(orth<=64*n*eps, "ev(true) orthogonality |V^H*V-I|", orth, 64*n*eps);

#ifdef WITH_REFERENCE_LAPACK
    for (int wantV=0; wantV<=1; ++wantV) {
        ZGeMatrix C = A;
        DVector      w(n), rw(std::max(1, 3*int(n)-2));
        int n_ = n, lWork = -1, info_;
        Z   wsq;
        zheev_(wantV ? "V" : "N", "U", &n_, C.data(), &n_, w.data(),
               &wsq, &lWork, rw.data(), &info_);
        lWork = int(real(wsq));
        ZVector wk(lWork);
        zheev_(wantV ? "V" : "N", "U", &n_, C.data(), &n_, w.data(),
               wk.data(), &lWork, rw.data(), &info_);
        check(info_==0, "zheev info", info_, 0);
        check(maxAbsDiff(wantV ? wV : wN, w)<=bound,
              wantV ? "ev(true) vs. reference zheev('V')"
                    : "ev(false) vs. reference zheev('N')",
              maxAbsDiff(wantV ? wV : wN, w), bound);
    }
#endif
}

int
main()
{
    // 2x2 and 3x3 cases from the original bug report
    testToeplitz(2, 1, 1);
    testToeplitz(3, 2, 1);

    for (Index n=1; n<=40; ++n) {
        testToeplitz(n, 2, -1);
        testToeplitz(n, -3, 0.5);
        testToeplitz(n, 0, 1e-3);
        testToeplitz(n, 1e8, 3e7);
    }
    for (T D=1e10; D<=1e300; D*=1e10) {
        testWideRange(D);
    }
    for (Index n=1; n<=25; ++n) {
        for (int seed=1; seed<=4; ++seed) {
            testHermitian(n, seed*1000+n);
        }
    }

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "sterf: all checks passed" << endl;
    return 0;
}
