#ifndef FLENS_LAPACK_INTERFACE_INLCUDE_CONFIG_H
#define FLENS_LAPACK_INTERFACE_INLCUDE_CONFIG_H 1

#include <complex>
#include <iomanip>
#include <iostream>

#ifdef LAPACK_DECL
#   undef   LAPACK_DECL
#endif
#define  LAPACK_DECL(x)     x##_

//
//  define typedefs for FLENS matrix and vector types
//
#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif
#include <flens/flens.cxx>

#ifndef INTEGER
#define INTEGER             int
#endif // INTEGER

#define DOUBLE              double
#define DOUBLE_COMPLEX      double
#define CXX_DOUBLE_COMPLEX  std::complex<double>


extern "C" {

void
LAPACK_DECL(xerbla)(const char *SRNAME, const int *INFO, int SRNAME_LEN);

} // extern "C"

#define LAPACK_ERROR(name, info)   LAPACK_DECL(xerbla)(name, info, strlen(name))

namespace flens {

template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Transpose>::value,
                    Transpose>::Type
convertTo(const char c)
{
    if (c=='N') {
        return NoTrans;
    } else if (c=='T') {
        return Trans;
    } else if (c=='C') {
        return ConjTrans;
    } else if (c=='R') {
        return Conj;
    }
    std::cerr << "c = " << c << std::endl;
    ASSERT(0);
    return NoTrans;
}


// matrix/vector types with DOUBLE
typedef FullStorageView<DOUBLE, cxxblas::ColMajor>       DFSView;
typedef ConstFullStorageView<DOUBLE, cxxblas::ColMajor>  DConstFSView;

typedef GeMatrix<DFSView>                                DGeMatrixView;
typedef GeMatrix<DConstFSView>                           DConstGeMatrixView;

typedef TrMatrix<DFSView>                                DTrMatrixView;
typedef TrMatrix<DConstFSView>                           DConstTrMatrixView;

typedef SyMatrix<DFSView>                                DSyMatrixView;
typedef SyMatrix<DConstFSView>                           DConstSyMatrixView;

typedef ArrayView<DOUBLE>                                DArrayView;
typedef DenseVector<DArrayView>                          DDenseVectorView;

typedef ConstArrayView<DOUBLE>                           DConstArrayView;
typedef DenseVector<DConstArrayView>                     DConstDenseVectorView;

// matrix/vector types with CXX_DOUBLE_COMPLEX
typedef FullStorageView<CXX_DOUBLE_COMPLEX,
                        cxxblas::ColMajor>               ZFSView;
typedef ConstFullStorageView<CXX_DOUBLE_COMPLEX,
                             cxxblas::ColMajor>          ZConstFSView;

typedef GeMatrix<ZFSView>                                ZGeMatrixView;
typedef GeMatrix<ZConstFSView>                           ZConstGeMatrixView;

typedef SyMatrix<ZFSView>                                ZSyMatrixView;
typedef SyMatrix<ZConstFSView>                           ZConstSyMatrixView;

typedef HeMatrix<ZFSView>                                ZHeMatrixView;
typedef HeMatrix<ZConstFSView>                           ZConstHeMatrixView;

typedef TrMatrix<ZFSView>                                ZTrMatrixView;
typedef TrMatrix<ZConstFSView>                           ZConstTrMatrixView;

typedef ArrayView<CXX_DOUBLE_COMPLEX>                    ZArrayView;
typedef DenseVector<ZArrayView>                          ZDenseVectorView;

typedef ConstArrayView<CXX_DOUBLE_COMPLEX>               ZConstArrayView;
typedef DenseVector<ZConstArrayView>                     ZConstDenseVectorView;

// matrix/vector types with bool
// TODO: using "bool" should be possible too
typedef Array<INTEGER>                                   BArray;
typedef ArrayView<INTEGER>                               BArrayView;
typedef DenseVector<BArray>                              BDenseVector;
typedef DenseVector<BArrayView>                          BDenseVectorView;

// matrix/vector types with INTEGER
typedef Array<INTEGER>                                   IArray;
typedef ArrayView<INTEGER>                               IArrayView;
typedef ConstArrayView<INTEGER>                          IConstArrayView;
typedef DenseVector<IArray>                              IDenseVector;
typedef DenseVector<IArrayView>                          IDenseVectorView;
typedef DenseVector<IConstArrayView>                     IConstDenseVectorView;


} // namespace flens

#endif // FLENS_LAPACK_INTERFACE_INCLUDE_CONFIG_H
