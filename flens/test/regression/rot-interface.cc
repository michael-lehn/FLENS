//
//  Regression test for the Fortran interface of the plane rotations in
//  blas_FLENS.a (flens/blas/interface/blas/rot.cc): srot_, drot_ (BLAS) and
//  crot_, zrot_ (LAPACK auxiliary routines CROT, ZROT).
//
//  CROT and ZROT take a real cosine C and a complex sine S.  crot_ and zrot_
//  of blas_FLENS.a declared C as complex: they read 8 (16) bytes where the
//  caller passes a real number of 4 (8) bytes, and with a complex C the
//  generic rotation y <- C*y - S*x was used instead of y <- C*y - conj(S)*x.
//
//  The routines compiled from the FLENS sources included below are renamed
//  to flens_srot_, ... and compared with srot_, ... of the BLAS/LAPACK
//  library (build with REFLAPACK=...).
//

#define BLAS(x)   flens_##x##_
#include <flens/blas/interface/blas/rot.cc>

#include <cxxstd/algorithm.h>
#include <cxxstd/cmath.h>
#include <cxxstd/complex.h>
#include <cxxstd/iostream.h>
#include <vector>

using namespace std;

extern "C" {

void srot_(const int *N, float *X, const int *INCX, float *Y, const int *INCY,
           const float *C, const float *S);

void drot_(const int *N, double *X, const int *INCX, double *Y,
           const int *INCY, const double *C, const double *S);

void crot_(const int *N, complex<float> *X, const int *INCX,
           complex<float> *Y, const int *INCY, const float *C,
           const complex<float> *S);

void zrot_(const int *N, complex<double> *X, const int *INCX,
           complex<double> *Y, const int *INCY, const double *C,
           const complex<double> *S);

} // extern "C"

static int failures = 0;

template <typename T>
T
randomValue(T *)
{
    return T(rand()%2001-1000)/T(250);
}

template <typename T>
complex<T>
randomValue(complex<T> *)
{
    return complex<T>(randomValue((T *)0), randomValue((T *)0));
}

//
//  Interface signatures as in the reference BLAS/LAPACK
//
template <typename TV, typename TC, typename TS>
struct Rot
{
    typedef void (*Function)(const int *, TV *, const int *, TV *, const int *,
                             const TC *, const TS *);
};

//
//  Apply the FLENS and the reference routine to the same data and compare.
//
template <typename TV, typename TC, typename TS>
void
compare(const char *name,
        typename Rot<TV,TC,TS>::Function flensRot,
        typename Rot<TV,TC,TS>::Function refRot,
        int n, int incX, int incY, TC c, TS s)
{
    typedef decltype(abs(TV(0)))  T;

    const T eps = numeric_limits<T>::epsilon();

    const int lenX = 1 + std::max(n-1, 0)*abs(incX);
    const int lenY = 1 + std::max(n-1, 0)*abs(incY);

//  one guard element on each side
    vector<TV> x(lenX+2), y(lenY+2);
    for (auto &v : x) v = randomValue((TV *)0);
    for (auto &v : y) v = randomValue((TV *)0);

    vector<TV> xRef = x, yRef = y;

    flensRot(&n, x.data()+1, &incX, y.data()+1, &incY, &c, &s);
    refRot(&n, xRef.data()+1, &incX, yRef.data()+1, &incY, &c, &s);

    const T scale = abs(c) + abs(s);
    T maxErr = 0;
    bool guardsOk = (x.front()==xRef.front()) && (x.back()==xRef.back())
                 && (y.front()==yRef.front()) && (y.back()==yRef.back());
    for (int i=0; i<lenX+2; ++i) {
        maxErr = std::max(maxErr, abs(x[i]-xRef[i]) / std::max(T(1), abs(xRef[i])));
    }
    for (int i=0; i<lenY+2; ++i) {
        maxErr = std::max(maxErr, abs(y[i]-yRef[i]) / std::max(T(1), abs(yRef[i])));
    }
    if (!guardsOk || maxErr>8*eps*std::max(T(1), scale)) {
        ++failures;
        cerr << "FAILED: " << name << " n=" << n << " incX=" << incX
             << " incY=" << incY << " c=" << c << " s=" << s
             << ": max rel. error " << maxErr
             << (guardsOk ? "" : " (guard element changed)") << endl;
    }
}

template <typename TV, typename TC, typename TS>
void
testRot(const char *name,
        typename Rot<TV,TC,TS>::Function flensRot,
        typename Rot<TV,TC,TS>::Function refRot,
        const vector<TS> &sines)
{
//  N<0 is not tested: like all level 1 routines of blas_FLENS.a these
//  routines abort for N<0, whereas the reference returns immediately.
    const int dims[] = { 0, 1, 2, 5, 9 };
    const int incs[] = { 1, 2, 3, -1, -2 };

    for (int n : dims) {
        for (int incX : incs) {
            for (int incY : incs) {
                for (const TS &s : sines) {
                    compare<TV,TC,TS>(name, flensRot, refRot,
                                      n, incX, incY, TC(0.6), s);
                }
            }
        }
    }
}

//
//  The C++ interface blas::rot is not changed by the fix of the Fortran
//  interface: with real c and complex s it is ZROT, with complex c and s it
//  is the generic rotation x <- c*x + s*y, y <- c*y - s*x.
//
static void
testCxxInterface()
{
    using namespace flens;

    typedef complex<double>  Z;

    const Z x0[3] = { Z(1,2), Z(3,-1), Z(-2,0.5) };
    const Z y0[3] = { Z(0.5,0), Z(-2,1), Z(4,-3) };

    ZDenseVector x(3), y(3);
    for (int i=1; i<=3; ++i) {
        x(i) = x0[i-1];
        y(i) = y0[i-1];
    }

    const double c = 0.6;
    const Z      s(0.48, 0.64);

    blas::rot(x, y, c, s);
    double err = 0;
    for (int i=1; i<=3; ++i) {
        err = std::max(err, abs(x(i) - (c*x0[i-1] + s*y0[i-1])));
        err = std::max(err, abs(y(i) - (c*y0[i-1] - conj(s)*x0[i-1])));
    }

    const Z cz(0.6, 0.1);
    for (int i=1; i<=3; ++i) {
        x(i) = x0[i-1];
        y(i) = y0[i-1];
    }
    blas::rot(x, y, cz, s);
    for (int i=1; i<=3; ++i) {
        err = std::max(err, abs(x(i) - (cz*x0[i-1] + s*y0[i-1])));
        err = std::max(err, abs(y(i) - (cz*y0[i-1] - s*x0[i-1])));
    }

    if (err>1e-14) {
        ++failures;
        cerr << "FAILED: blas::rot (C++ interface): error " << err << endl;
    }
}

int
main()
{
    srand(4711);

    testRot<float,float,float>("srot_", flens_srot_, srot_,
                               { 0.8f, -0.8f, 0.0f, 2.5f });
    testRot<double,double,double>("drot_", flens_drot_, drot_,
                                  { 0.8, -0.8, 0.0, 2.5 });

    typedef complex<float>   C;
    typedef complex<double>  Z;

    testRot<C,float,C>("crot_", flens_crot_, crot_,
                       { C(0.8f,0), C(0,0.8f), C(0.48f,-0.64f), C(1.5f,2) });
    testRot<Z,double,Z>("zrot_", flens_zrot_, zrot_,
                        { Z(0.8,0), Z(0,0.8), Z(0.48,-0.64), Z(1.5,2) });

    testCxxInterface();

    if (failures) {
        cerr << failures << " failure(s)" << endl;
        return 1;
    }
    cout << "rot-interface: all tests passed" << endl;
    return 0;
}
