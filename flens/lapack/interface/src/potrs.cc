#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dpotrs --------------------------------------------------------------------
void
LAPACK_DECL(dpotrs)(const char       *UPLO,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
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
        LAPACK_ERROR("DPOTRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo         upLo = StorageUpLo(*UPLO);
    DConstFSView        AFS  = DConstFSView(*N, *N, A, *LDA);
    DConstSyMatrixView  _A   = DConstSyMatrixView(AFS, upLo);
    DGeMatrixView       _B   = DFSView(*N, *NRHS, B, *LDB);

    potrs(_A, _B);
}

//-- zpotrs --------------------------------------------------------------------
void
LAPACK_DECL(zpotrs)(const char               *UPLO,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
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
        LAPACK_ERROR("ZPOTRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto zA = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto zB       = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);

    StorageUpLo         upLo = StorageUpLo(*UPLO);
    ZConstFSView        AFS  = ZConstFSView(*N, *N, zA, *LDA);
    ZConstHeMatrixView  _A   = ZConstHeMatrixView(AFS, upLo);
    ZGeMatrixView       _B   = ZFSView(*N, *NRHS, zB, *LDB);

    potrs(_A, _B);
}

} // extern "C"

} } // namespace lapack, flens
