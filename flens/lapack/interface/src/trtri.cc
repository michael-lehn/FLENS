#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dtrtri --------------------------------------------------------------------
void
LAPACK_DECL(dtrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool upper  = (*UPLO=='U');
    const bool noUnit = (*DIAG=='N');

    if (!upper && *UPLO!='L') {
        *INFO = 1;
    } else if (!noUnit && *DIAG!='U') {
        *INFO = 2;
    } else if (*N<0) {
        *INFO = 3;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 5;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DTRTRI", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    StorageUpLo    upLo = cxxblas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Diag           diag = cxxblas::getCxxBlasEnum<Diag>(*DIAG);

    DTrMatrixView  A_(DFSView(*N, *N, A, *LDA), upLo, diag);

//
//  Call FLENS implementation
//
    tri(A_);
}

//-- ztrtri --------------------------------------------------------------------
void
LAPACK_DECL(ztrtri)(const char       *UPLO,
                    const char       *DIAG,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *INFO)
{
    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool upper  = (*UPLO=='U');
    const bool noUnit = (*DIAG=='N');

    if (!upper && *UPLO!='L') {
        *INFO = 1;
    } else if (!noUnit && *DIAG!='U') {
        *INFO = 2;
    } else if (*N<0) {
        *INFO = 3;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 5;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZTRTRI", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    StorageUpLo    upLo = cxxblas::getCxxBlasEnum<StorageUpLo>(*UPLO);
    Diag           diag = cxxblas::getCxxBlasEnum<Diag>(*DIAG);

    auto zA = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    ZTrMatrixView  A_(ZFSView(*N, *N, zA, *LDA), upLo, diag);

//
//  Call FLENS implementation
//
    tri(A_);
}

} // extern "C"

} } // namespace lapack, flens
