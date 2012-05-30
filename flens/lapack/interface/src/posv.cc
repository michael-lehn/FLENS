//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dposv ---------------------------------------------------------------------
void
LAPACK_DECL(dposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("dposv");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*NRHS<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -5;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DPOSV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo    upLo = StorageUpLo(*UPLO);
    DFSView        AFS  = DFSView(*N, *N, A, *LDA);
    DSyMatrixView  _A   = DSyMatrixView(AFS, upLo);
    DGeMatrixView  _B   = DFSView(*N, *NRHS, B, *LDB);

    posv(_A, _B);
}

//-- zposv ---------------------------------------------------------------------
void
LAPACK_DECL(zposv)(const char           *UPLO,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("zposv");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*NRHS<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -5;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZPOSV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zB     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);

    StorageUpLo    upLo = StorageUpLo(*UPLO);
    ZFSView        AFS  = ZFSView(*N, *N, zA, *LDA);
    ZHeMatrixView  _A   = ZHeMatrixView(AFS, upLo);
    ZGeMatrixView  _B   = ZFSView(*N, *NRHS, zB, *LDB);

    posv(_A, _B);
}

} // extern "C"

} } // namespace lapack, flens
