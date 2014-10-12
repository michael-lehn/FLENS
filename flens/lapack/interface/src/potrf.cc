#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dpotrf --------------------------------------------------------------------
void
LAPACK_DECL(dpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
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
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DPOTRF", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo    upLo = StorageUpLo(*UPLO);
    DFSView        AFS  = DFSView(*N, *N, A, *LDA);
    DSyMatrixView  A_   = DSyMatrixView(AFS, upLo);

    *INFO = potrf(A_);
}

//-- zpotrf --------------------------------------------------------------------
void
LAPACK_DECL(zpotrf)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
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
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZPOTRF", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    auto zA = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);

    StorageUpLo    upLo = StorageUpLo(*UPLO);
    ZFSView        AFS  = ZFSView(*N, *N, zA, *LDA);
    ZHeMatrixView  A_   = ZHeMatrixView(AFS, upLo);

    *INFO = potrf(A_);
}

} // extern "C"

} } // namespace lapack, flens
