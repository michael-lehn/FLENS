#define STR(x)      #x
#define STRING(x)   STR(x)

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
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dtrtrs");

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
    StorageUpLo         upLo   = cxxf77blas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Transpose           trans  = cxxf77blas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxf77blas::getCxxBlasEnum<Diag>(*DIAG);
    DConstTrMatrixView  _A(DConstFSView(*N, *N, *LDA, A), upLo, diag);
    DGeMatrixView       _B     = DFSView(*N, *NRHS, *LDB, B);

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
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: ztrtrs");

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

    StorageUpLo         upLo   = cxxf77blas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Transpose           trans  = cxxf77blas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxf77blas::getCxxBlasEnum<Diag>(*DIAG);
    ZConstTrMatrixView  _A(ZConstFSView(*N, *N, *LDA, zA), upLo, diag);
    ZGeMatrixView       _B     = ZFSView(*N, *NRHS, *LDB, zB);

    trs(trans, _A, _B);
}

} // extern "C"

} } // namespace lapack, flens
