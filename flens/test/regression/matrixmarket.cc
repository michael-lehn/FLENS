//
//  Test for loadMatrixMarket.
//
//  When assertions are enabled, loadMatrixMarket converts the header line
//  to lower case and checks it.  This conversion used std::locale(""), the
//  user's preferred locale from the environment, which throws if the locale
//  named in the environment is not available.  With libstdc++ on macOS this
//  happens for any locale except "C"/"POSIX" (e.g. LANG=en_US.UTF-8), so
//  loading a file terminated the program.  The test therefore sets LC_ALL
//  to an invalid locale name before loading.
//
//  Files with a mixed case header ("%%MatrixMarket ...") are loaded into
//  GeMatrix and SpMatrix (real and complex).  The GeMatrix files are written
//  with saveMatrixMarket.
//
//  Known problems not covered here (they are independent of the C++
//  standard and exist in the original code):
//    - loadMatrixMarket for SyMatrix and HeMatrix calls
//      A.resize(numRows, numCols), i.e. resize(dim, firstIndex), and then
//      fails with an index assertion.  These loaders are only instantiated
//      (see compileOnly) to make sure they compile.
//    - saveMatrixMarket for SpMatrix/HpMatrix and loadMatrixMarket for
//      HpMatrix do not compile (they use firstIndex()/lastIndex(), which
//      the packed matrix types do not provide).  The SpMatrix files are
//      therefore written directly.
//
//  The test writes temporary files into the current directory.
//

#include <cxxstd/complex.h>
#include <cxxstd/cstdio.h>
#include <cxxstd/fstream.h>
#include <cxxstd/iostream.h>
#include <flens/flens.cxx>

#include <stdlib.h>      // setenv (POSIX)

using namespace flens;
using namespace std;

typedef complex<double>  Z;

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
//  exactly representable test values
//
template <typename T>
static T
value(int i, int j);

template <>
double
value<double>(int i, int j)
{
    return (8*i - 3*j)/8.0;
}

template <>
Z
value<Z>(int i, int j)
{
    return Z(value<double>(i, j), (i+j)/4.0);
}

template <typename T>
static void
testGe(const char *file)
{
    GeMatrix<FullStorage<T> > A(4, 3), B;
    for (int i=1; i<=4; ++i) {
        for (int j=1; j<=3; ++j) {
            A(i,j) = value<T>(i, j);
        }
    }
    check(saveMatrixMarket(file, A), "saveMatrixMarket(GeMatrix)");
    check(loadMatrixMarket(file, B), "loadMatrixMarket(GeMatrix)");
    check(B.numRows()==4 && B.numCols()==3, "GeMatrix: dimensions");
    for (int i=1; i<=4; ++i) {
        for (int j=1; j<=3; ++j) {
            check(B(i,j)==A(i,j), "GeMatrix: element");
        }
    }
    remove(file);
}

//
//  Instantiate (but do not call) the loaders for SyMatrix and HeMatrix.
//
void
compileOnly(const char *file)
{
    SyMatrix<FullStorage<double> >  S(1, Upper);
    SyMatrix<FullStorage<Z> >       SZ(1, Upper);
    HeMatrix<FullStorage<Z> >       H(1, Upper);
    loadMatrixMarket(file, S);
    loadMatrixMarket(file, SZ);
    loadMatrixMarket(file, H);
}

static void
writeEntry(ofstream &ofs, double a)
{
    ofs << a << endl;
}

static void
writeEntry(ofstream &ofs, const Z &a)
{
    ofs << real(a) << " " << imag(a) << endl;
}

template <typename T>
static void
testSp(const char *file, const char *field)
{
    // upper triangle row by row, as expected by loadMatrixMarket(SpMatrix)
    {
        ofstream ofs(file);
        ofs << "%%MatrixMarket Matrix Array " << field << " Symmetric" << endl;
        ofs << "% written by the FLENS regression test" << endl;
        ofs << "4 4" << endl;
        for (int i=1; i<=4; ++i) {
            for (int j=i; j<=4; ++j) {
                writeEntry(ofs, value<T>(i, j));
            }
        }
    }
    SpMatrix<PackedStorage<T> > B(1);
    check(loadMatrixMarket(file, B), "loadMatrixMarket(SpMatrix)");
    check(B.dim()==4, "SpMatrix: dimension");
    for (int i=1; i<=4; ++i) {
        for (int j=i; j<=4; ++j) {
            check(B(i,j)==value<T>(i, j), "SpMatrix: element");
        }
    }
    remove(file);
}

int
main()
{
    // header parsing must not depend on the locale of the environment
    setenv("LC_ALL", "flens-test-invalid-locale", 1);

    testGe<double>("matrixmarket-ge-real.mtx");
    testGe<Z>("matrixmarket-ge-complex.mtx");
    testSp<double>("matrixmarket-sp-real.mtx", "Real");
    testSp<Z>("matrixmarket-sp-complex.mtx", "Complex");

    if (failures) {
        cerr << failures << " check(s) FAILED" << endl;
        return 1;
    }
    cout << "matrixmarket: all checks passed" << endl;
    return 0;
}
