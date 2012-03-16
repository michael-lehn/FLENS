#include <complex>
#include <iostream>

///
/// Define `USE_CXXLAPACK` before you include the FLENS headers
///
#define USE_CXXLAPACK
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

///
/// Elements are now defined to be *complex*
///
typedef complex<double>   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;
    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 3;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    A = T(1, 0), T( 1,-1), T(  2,20),
        T(0, 1), T( 1, 2), T(-10, 8),
        T(0,-1), T(-1, 1), T( 40, 6);

    b = T( 1, 0),
        T(-1, 1),
        T(-2,-1);

    cerr << "A = " << A << endl;
    cerr << "b = " << b << endl;

///
/// We use __cxxlapack::getrf__ as we did in the previous example ...
///
/// :links: __cxxlapack::getrf__ -> file:cxxlapack/interface/getrf.h
///
    IndexType info = cxxlapack::getrf(A.numRows(),
                                      A.numCols(),
                                      A.data(),
                                      A.leadingDimension(),
                                      piv.data());

    if (info==0) {
///
///     ... and also __cxxlapack::getrs__.
///
///     :links: __cxxlapack::getrs__ -> file:cxxlapack/interface/getrs.h
///
        cxxlapack::getrs(NoTrans,
                         A.numRows(),
                         IndexType(1),
                         A.data(),
                         A.leadingDimension(),
                         piv.data(),
                         b.data(),
                         b.length());

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }

}

