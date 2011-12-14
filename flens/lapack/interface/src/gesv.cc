//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgesv ---------------------------------------------------------------------
void
LAPACK_DECL(dgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("dgesv");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*N<0) {
        *INFO = -1;
    } else if (*NRHS<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGESV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    DGeMatrixView     _A     = DFSView(*N, *N, A, *LDA);
    IDenseVectorView  _IPIV  = IArrayView(*N, IPIV, 1);
    DGeMatrixView     _B     = DFSView(*N, *NRHS, B, *LDB);

    sv(_A, _IPIV, _B);
}

} // extern "C"

} } // namespace lapack, flens
