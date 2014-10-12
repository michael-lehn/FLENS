#define STR(x)      #x
#define STRING(x)   STR(x)

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
    DGeMatrixView     A_     = DFSView(*N, *N, A, *LDA);
    IDenseVectorView  IPIV_  = IArrayView(*N, IPIV, 1);
    DGeMatrixView     B_     = DFSView(*N, *NRHS, B, *LDB);

    sv(A_, IPIV_, B_);
}

//-- zgesv ---------------------------------------------------------------------
void
LAPACK_DECL(zgesv)(const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   INTEGER              *IPIV,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   INTEGER              *INFO)
{
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
        LAPACK_ERROR("ZGESV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zB     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);

    ZGeMatrixView     A_     = ZFSView(*N, *N, zA, *LDA);
    IDenseVectorView  IPIV_  = IArrayView(*N, IPIV, 1);
    ZGeMatrixView     B_     = ZFSView(*N, *NRHS, zB, *LDB);

    sv(A_, IPIV_, B_);
}

} // extern "C"

} } // namespace lapack, flens
