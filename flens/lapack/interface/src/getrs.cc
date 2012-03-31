//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgetrs --------------------------------------------------------------------
void
LAPACK_DECL(dgetrs)(const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgetrs");
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
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -8;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGETRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    Transpose              trans  = convertTo<Transpose>(*TRANS);
    DConstGeMatrixView     _A     = DConstFSView(*N, *N, A, *LDA);
    IConstDenseVectorView  _IPIV  = IConstArrayView(*N, IPIV, 1);
    DGeMatrixView          _B     = DFSView(*N, *NRHS, B, *LDB);

    trs(trans, _A, _IPIV, _B);
}

} // extern "C"

} } // namespace lapack, flens
