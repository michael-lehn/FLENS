//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dtrtrs --------------------------------------------------------------------
void
LAPACK_DECL(dtrtrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("--> dtrtrs");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;

    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*TRANS!='N' && *TRANS!='C' && *TRANS!='T') {
        *INFO = -2;
    } else if (*DIAG!='N' && *DIAG!='U') {
        *INFO = -3;
    } else if (*N<0) {
        *INFO = -4;
    } else if (*NRHS<0) {
        *INFO = -5;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -9;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DTRTRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo         upLo   = cxxblas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Transpose           trans  = cxxblas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxblas::getCxxBlasEnum<Diag>(*DIAG);
    DConstTrMatrixView  _A(DConstFSView(*N, *N, A, *LDA), upLo, diag);
    DGeMatrixView       _B     = DFSView(*N, *NRHS, B, *LDB);

    trs(trans, _A, _B);
}

//-- ztrtrs --------------------------------------------------------------------
void
LAPACK_DECL(ztrtrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
    DEBUG_FLENS_LAPACK("--> ztrtrs");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;

    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*TRANS!='N' && *TRANS!='C' && *TRANS!='T') {
        *INFO = -2;
    } else if (*DIAG!='N' && *DIAG!='U') {
        *INFO = -3;
    } else if (*N<0) {
        *INFO = -4;
    } else if (*NRHS<0) {
        *INFO = -5;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -9;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZTRTRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto zA = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto zB       = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);

    StorageUpLo         upLo   = cxxblas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Transpose           trans  = cxxblas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxblas::getCxxBlasEnum<Diag>(*DIAG);
    ZConstTrMatrixView  _A(ZConstFSView(*N, *N, zA, *LDA), upLo, diag);
    ZGeMatrixView       _B     = ZFSView(*N, *NRHS, zB, *LDB);

    trs(trans, _A, _B);
}

} // extern "C"

} } // namespace lapack, flens
