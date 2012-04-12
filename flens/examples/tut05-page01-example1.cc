#include <iostream>

///
/// Define `USE_CXXLAPACK` before you include the FLENS headers
///
#define USE_CXXLAPACK
#include <flens/flens.cxx>

using namespace std;
using namespace flens;

typedef double   T;

int
main()
{
    typedef GeMatrix<FullStorage<T> >           Matrix;
    typedef DenseVector<Array<T> >              Vector;

    typedef Matrix::IndexType                   IndexType;
    typedef DenseVector<Array<IndexType> >      IndexVector;

    const IndexType n = 4;

    Matrix         A(n,n);
    Vector         b(n);
    IndexVector    piv(n);

    A =  2,   3,  -1,   0,
        -6,  -5,   0,   2,
         2,  -5,   6,  -6,
         4,   6,   2,  -3;

    b = 20,
       -33,
       -43,
        49;

    cerr << "A = " << A << endl;
    cerr << "b = " << b << endl;

///
/// The CXXLAPACK functions require the usual information describing the
/// matrices, i.e. dimensions, leading dimension and a data pointer.  All
/// these can be retrieved from the FLENS matrix type.
///
    IndexType info = cxxlapack::getrf(A.numRows(),
                                      A.numCols(),
                                      A.data(),
                                      A.leadingDimension(),
                                      piv.data());

    if (info==0) {
///
///     Note that __cxxlapack::getrs__ has a template parameter for the used
///     index type.  That's why we call it with `IndexType(1)`.  This makes
///     sure that values `A.numRows()`, `IndexType(1)`, `A.leadingDimension`
///     and `b.stride()` all have the same type.
///
///     :links: __cxxlapack::getrs__ -> file:cxxlapack/interface/getrs.h
///
        cxxlapack::getrs(lapack::getF77Char(NoTrans),
                         A.numRows(),
                         IndexType(1),
                         A.data(),
                         A.leadingDimension(),
                         piv.data(),
                         b.data(),
                         b.length());

///
///     As an alternative we could specify the template parameter explictly:
///
/*
        cxxlapack::getrs<IndexType>(lapack::getF77Char(NoTrans),
                                    A.numRows(),
                                    1,
                                    A.data(),
                                    A.leadingDimension(),
                                    piv.data(),
                                    b.data(),
                                    b.stride());
*/

        cerr << "x = " << b << endl;
    } else {
        cerr << "A is numerically singular." << endl;
    }
}

