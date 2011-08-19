#ifndef FLENS_LAPACK_INTERFACE_INTERFACE_H
#define FLENS_LAPACK_INTERFACE_INTERFACE_H 1


//#define CXXBLAS_DEBUG_OUT(msg)      std::cerr << msg << std::endl;



#include <cassert>
#include <cmath>
#include <iostream>

using std::abs;

//------------------------------------------------------------------------------
//
//  Use ilaenv from lapack test suite.
//  Note: These functions have to be defined before including flens.cxx
//

extern "C" {

#define ILAENV          ilaenv_

int
ILAENV(int *SPEC, const char *NAME, const char *OPTS,
       int *N1, int *N2, int *N3, int *N4);

} // extern "C"

#define FLENS_LAPACK_AUX_ILAENV_H 1
#define FLENS_LAPACK_AUX_ILAENV_TCC 1

#include <flens/aux/issame.h>
#include <string>

namespace flens { namespace lapack {

using std::string;

template <typename T>
int
ilaenv(int spec, const char *_name, const char *_opts,
       int n1=-1, int n2=-1, int n3=-1, int n4=-1)
{
    /*
    std::cerr << "spec = " << spec
              << ", name = " << _name
              << ", opts = " << _opts
              << ", n1 = " << n1
              << ", n2 = " << n2
              << ", n3 = " << n3
              << ", n4 = " << n4
              << std::endl;
    */

    string opts(_opts);
    string name("D");
    name = name + string(_name);

    int result = ILAENV(&spec, name.c_str(), opts.c_str(), &n1, &n2, &n3, &n4);

    /*
    std::cerr << " => result = " << result << std::endl;
    */

    return result;
}

} } // namespace lapack, flens

//
//------------------------------------------------------------------------------
//
//  Define element and index types
//

#ifndef INT
#   define INT                      int
#endif

#ifndef FLENS_DEFAULT_INDEXTYPE
#   define FLENS_DEFAULT_INDEXTYPE  int
#endif

#include <cstring>
#include <flens/flens.cxx>

#ifdef SINGLE
#   define FLOAT        float
#   define CXXFLOAT     float
#elif DOUBLE
#   define FLOAT        double
#   define CXXFLOAT     double
#elif COMPLEX_SINGLE
#   define FLOAT        float
#   define CXXFLOAT     std::complex<float>
#elif COMPLEX_DOUBLE
#   define FLOAT        double
#   define CXXFLOAT     std::complex<double>
#endif

//
//------------------------------------------------------------------------------
//
//  Import common types and enums.
//

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;

using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;

using cxxblas::Side;
using cxxblas::Left;
using cxxblas::Right;

using cxxblas::Transpose;
using cxxblas::NoTrans;
using cxxblas::Trans;
using cxxblas::Conj;
using cxxblas::ConjTrans;

using cxxblas::Diag;
using cxxblas::Unit;
using cxxblas::NonUnit;

//
//------------------------------------------------------------------------------
//
//  Define common matrix and vector types
//

namespace flens {

typedef IndexOptions<INT>           IdxOpt;
typedef std::allocator<FLOAT>       Alloc;

typedef FullStorageView<CXXFLOAT, cxxblas::ColMajor, IdxOpt, Alloc>  FSV;
typedef FullStorage<CXXFLOAT, cxxblas::ColMajor>                     FS;

typedef Array<CXXFLOAT, IdxOpt, Alloc>                               AR;
typedef ArrayView<CXXFLOAT, IdxOpt, Alloc>                           AV;
typedef ArrayView<INT, IdxOpt, Alloc>                                IAV;
typedef Array<INT>                                                   IA;

} // namespace flens

//
//------------------------------------------------------------------------------
//
//  Functions for comparing results
//

namespace flens {

template <typename MA, typename MB>
bool
isDifferent(const GeMatrix<MA> &A, const GeMatrix<MB> &B)
{
    typedef typename GeMatrix<MA>::IndexType IndexType;

    if (A.numRows()!=B.numRows()) {
        std::cerr << "A.numRows() = " << A.numRows()
                  << ", B.numRows() = " << B.numRows()
                  << std::endl;
        return true;
    }
    if (A.numCols()!=B.numCols()) {
        std::cerr << "A.numCols() = " << A.numCols()
                  << ", B.numCols() = " << B.numCols()
                  << std::endl;
        return true;
    }
    if (A.firstRow()!=B.firstRow()) {
        std::cerr << "A.firstRow() = " << A.firstRow()
                  << ", B.firstRow() = " << B.firstRow()
                  << std::endl;
        return true;
    }
    if (A.firstCol()!=B.firstCol()) {
        std::cerr << "A.firstCol() = " << A.firstCol()
                  << ", B.firstCol() = " << B.firstCol()
                  << std::endl;
        return true;
    }
    if (A.lastRow()!=B.lastRow()) {
        std::cerr << "A.lastRow() = " << A.lastRow()
                  << ", B.lastRow() = " << B.lastRow()
                  << std::endl;
        return true;
    }
    if (A.lastCol()!=B.lastCol()) {
        std::cerr << "A.lastCol() = " << A.lastCol()
                  << ", B.lastCol() = " << B.lastCol()
                  << std::endl;
        return true;
    }

    for (IndexType i=A.firstRow(); i<=A.lastRow(); ++i) {
        for (IndexType j=A.firstCol(); j<=A.lastCol(); ++j) {
            if (A(i,j)!=B(i,j)) {
                std::cerr.precision(150);
                std::cerr << "A(" << i << ", " << j << ") = " << A(i,j)
                          << std::endl
                          << "B(" << i << ", " << j << ") = " << B(i,j)
                          << std::endl
                          << "A(" << i << ", " << j << ") - "
                          << "B(" << i << ", " << j << ") = "
                          << A(i,j) - B(i,j)
                          << std::endl;
                return true;
            }
        }
    }
    return false;
}


template <typename VX, typename VY>
bool
isDifferent(const DenseVector<VX> &x, const DenseVector<VY> &y)
{
    typedef typename DenseVector<VX>::IndexType IndexType;

    if (x.length()!=y.length()) {
        std::cerr << "x.length() = " << x.length()
                  << ", y.length() = " << y.length()
                  << std::endl;
        return true;
    }
    if (x.firstIndex()!=y.firstIndex()) {
        std::cerr << "x.firstIndex() = " << x.firstIndex()
                  << ", x.firstIndex() = " << x.firstIndex()
                  << std::endl;
        return true;
    }
    if (x.endIndex()!=y.endIndex()) {
        std::cerr << "x.endIndex() = " << x.endIndex()
                  << ", y.endIndex() = " << y.endIndex()
                  << std::endl;
        return true;
    }
    if (x.inc()!=y.inc()) {
        std::cerr << "x.inc() = " << x.inc()
                  << ", x.inc() = " << x.inc()
                  << std::endl;
        return true;
    }
    
    for (IndexType i=x.firstIndex(); i!=x.endIndex(); i+=x.inc()) {
        if (x(i)!=y(i)) {
            std::cerr.precision(150);
            std::cerr << "x(" << i << ") = " << x(i)
                      << std::endl
                      << "y(" << i << ") = " << y(i)
                      << std::endl
                      << "x(" << i << ") - "
                      << "y(" << i << ") = "
                      << x(i) - y(i)
                      << std::endl;
            return true;
        }
    }
    return false;
}



} // namespace flens

//
//------------------------------------------------------------------------------
//
//  Use xerbla from lapack test suite
//

extern "C" {

#define XERBLA          xerbla_

void
XERBLA(const char *SRNAME, int *INFO, int len);

} // extern "C"


#endif // FLENS_LAPACK_INTERFACE_INTERFACE_H
