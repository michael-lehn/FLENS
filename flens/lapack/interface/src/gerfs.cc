//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgerfs --------------------------------------------------------------------
void
LAPACK_DECL(dgerfs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *AF,
                    const INTEGER    *LDAF,
                    const INTEGER    *IPIV,
                    const DOUBLE     *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgerfs");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*TRANS!='N' && *TRANS!='T' && *TRANS!='C') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*NRHS<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -5;
    } else if (*LDAF<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -10;
    } else if (*LDX<std::max(INTEGER(1), *N)) {
        *INFO = -12;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGERFS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    Transpose              trans  = getFlensLapackEnum<Transpose>(*TRANS);
    DConstGeMatrixView     _A     = DConstFSView(*N, *N, A, *LDA);
    DConstGeMatrixView     _AF    = DConstFSView(*N, *N, AF, *LDAF);
    IConstDenseVectorView  _IPIV  = IConstArrayView(*N, IPIV, 1);
    DConstGeMatrixView     _B     = DConstFSView(*N, *NRHS, B, *LDB);
    DGeMatrixView          _X     = DFSView(*N, *NRHS, X, *LDX);
    DDenseVectorView       _FERR  = DArrayView(*NRHS, FERR, 1);
    DDenseVectorView       _BERR  = DArrayView(*NRHS, BERR, 1);
    DDenseVectorView       _WORK  = DArrayView(*N*3, WORK, 1);
    IDenseVectorView       _IWORK = IArrayView(*N, IWORK, 1);

    rfs(trans, _A, _AF, _IPIV, _B, _X, _FERR, _BERR, _WORK, _IWORK);
}

} // extern "C"

} } // namespace lapack, flens
