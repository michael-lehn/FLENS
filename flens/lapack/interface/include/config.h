#ifndef FLENS_LAPACK_INTERFACE_INLCUDE_CONFIG_H
#define FLENS_LAPACK_INTERFACE_INLCUDE_CONFIG_H 1

// #define CXXBLAS_DEBUG_OUT(x)    std::cerr << x << std::endl;


#include <complex>
#include <iomanip>
#include <iostream>

#ifndef INTEGER
#   define INTEGER          int
#endif

#define UNKNOWN             void

#define DOUBLE              double
#define FLOAT               float
#define LOGICAL             int
#define DOUBLE_COMPLEX      std::complex<double>

#include <flens/lapack/interface/include/cxxlapack.h>
#include <flens/lapack/interface/include/f77lapack.h>

#ifdef LAPACK_DECL
#   undef   LAPACK_DECL
#endif 
#define  LAPACK_DECL(x)     x##_

#ifdef LAPACK_IMPL
#   undef   LAPACK_IMPL
#endif
#define  LAPACK_IMPL(x)     x


#define DEBUG_LAPACK_CALL(msg)      std::cerr << "LAPACK_CALL: " \
                                              << msg << std::endl;

#ifdef DEBUG_LAPACK_STUB
#   define DEBUG_LAPACK_STUB(msg)      std::cerr << "LAPACK_STUB: " \
                                                 << msg << std::endl;
#else
#   define DEBUG_LAPACK_STUB(msg)
#endif


#define LAPACK_ERROR(name, info)    xerbla_(name, info, strlen(name));

//
//  define typedefs for FLENS matrix and vector types
//
#include <flens/flens.cxx>

namespace flens {

typedef FullStorageView<DOUBLE, cxxblas::ColMajor>  DFSView;
typedef GeMatrix<DFSView>                           DGeMatrixView;

typedef Array<bool>                                 BArray;
typedef ArrayView<bool>                             BArrayView;
typedef DenseVector<BArray>                         BDenseVector;
typedef DenseVector<BArrayView>                     BDenseVectorView;

typedef Array<INTEGER>                              IArray;
typedef ArrayView<INTEGER>                          IArrayView;
typedef DenseVector<IArray>                         IDenseVector;
typedef DenseVector<IArrayView>                     IDenseVectorView;

typedef ArrayView<DOUBLE>                           DArrayView;
typedef DenseVector<DArrayView>                     DDenseVectorView;

typedef ConstArrayView<DOUBLE>                      DConstArrayView;
typedef DenseVector<DConstArrayView>                DConstDenseVectorView;

} // namespace flens

#endif // FLENS_LAPACK_INTERFACE_INCLUDE_CONFIG_H
