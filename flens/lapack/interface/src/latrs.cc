#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dlatrs --------------------------------------------------------------------
void
LAPACK_DECL(dlatrs)(const char       *UPLO,
                    const char       *TRANS,
                    const char       *DIAG,
                    const char       *NORMIN,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *X,
                    DOUBLE           *SCALE,
                    DOUBLE           *CNORM,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dlatrs");

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
    } else if (*NORMIN!='Y' && *NORMIN!='N') {
        *INFO = -4;
    } else if (*N<0) {
        *INFO = -5;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DLATRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const bool          normIn = (*NORMIN=='Y') ? true : false;
    Transpose           trans  = cxxf77blas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxf77blas::getCxxBlasEnum<Diag>(*DIAG);
    StorageUpLo         upLo   = cxxf77blas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    DConstTrMatrixView  _A(DConstFSView(*N, *N, *LDA, A), upLo, diag);
    DDenseVectorView    _X     = DArrayView(*N, X, 1);
    DDenseVectorView    _CNORM = DArrayView(*N, CNORM, 1);

    latrs(trans, normIn, _A, _X, *SCALE, _CNORM);
}

//-- zlatrs --------------------------------------------------------------------
void
LAPACK_DECL(zlatrs)(const char               *UPLO,
                    const char               *TRANS,
                    const char               *DIAG,
                    const char               *NORMIN,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE_COMPLEX           *X,
                    DOUBLE                   *SCALE,
                    DOUBLE                   *CNORM,
                    INTEGER                  *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: zlatrs");

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
    } else if (*NORMIN!='Y' && *NORMIN!='N') {
        *INFO = -4;
    } else if (*N<0) {
        *INFO = -5;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZLATRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto zA   = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto zX         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(X);

    const bool          normIn = (*NORMIN=='Y') ? true : false;
    Transpose           trans  = cxxf77blas::getCxxBlasEnum<Transpose>(*TRANS);
    Diag                diag   = cxxf77blas::getCxxBlasEnum<Diag>(*DIAG);
    StorageUpLo         upLo   = cxxf77blas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    ZConstTrMatrixView  _A(ZConstFSView(*N, *N, *LDA, zA), upLo, diag);
    ZDenseVectorView    _X     = ZArrayView(*N, zX, 1);
    DDenseVectorView    _CNORM = DArrayView(*N, CNORM, 1);

    latrs(trans, normIn, _A, _X, *SCALE, _CNORM);
}

} // extern "C"

} } // namespace lapack, flens
