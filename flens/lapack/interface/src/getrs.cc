#define STR(x)      #x
#define STRING(x)   STR(x)

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
    DConstGeMatrixView     A_     = DConstFSView(*N, *N, A, *LDA);
    IConstDenseVectorView  IPIV_  = IConstArrayView(*N, IPIV, 1);
    DGeMatrixView          B_     = DFSView(*N, *NRHS, B, *LDB);

    trs(trans, A_, IPIV_, B_);
}

//-- zgetrs --------------------------------------------------------------------
void
LAPACK_DECL(zgetrs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const INTEGER            *IPIV,
                    DOUBLE_COMPLEX           *B,
                    const INTEGER            *LDB,
                    INTEGER                  *INFO)
{
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
        LAPACK_ERROR("ZGETRS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto zA = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto zB       = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);

    Transpose              trans  = convertTo<Transpose>(*TRANS);
    ZConstGeMatrixView     A_     = ZConstFSView(*N, *N, zA, *LDA);
    IConstDenseVectorView  IPIV_  = IConstArrayView(*N, IPIV, 1);
    ZGeMatrixView          B_     = ZFSView(*N, *NRHS, zB, *LDB);
    trs(trans, A_, IPIV_, B_);
}

} // extern "C"

} } // namespace lapack, flens
