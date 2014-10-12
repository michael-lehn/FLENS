#define STR(x)      #x
#define STRING(x)   STR(x)

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
    Transpose              trans  = convertTo<Transpose>(*TRANS);
    DConstGeMatrixView     A_     = DConstFSView(*N, *N, A, *LDA);
    DConstGeMatrixView     AF_    = DConstFSView(*N, *N, AF, *LDAF);
    IConstDenseVectorView  IPIV_  = IConstArrayView(*N, IPIV, 1);
    DConstGeMatrixView     B_     = DConstFSView(*N, *NRHS, B, *LDB);
    DGeMatrixView          X_     = DFSView(*N, *NRHS, X, *LDX);
    DDenseVectorView       FERR_  = DArrayView(*NRHS, FERR, 1);
    DDenseVectorView       BERR_  = DArrayView(*NRHS, BERR, 1);
    DDenseVectorView       WORK_  = DArrayView(*N*3, WORK, 1);
    IDenseVectorView       IWORK_ = IArrayView(*N, IWORK, 1);

    rfs(trans, A_, AF_, IPIV_, B_, X_, FERR_, BERR_, WORK_, IWORK_);
}

//-- zgerfs --------------------------------------------------------------------
void
LAPACK_DECL(zgerfs)(const char               *TRANS,
                    const INTEGER            *N,
                    const INTEGER            *NRHS,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *AF,
                    const INTEGER            *LDAF,
                    const INTEGER            *IPIV,
                    const DOUBLE_COMPLEX     *B,
                    const INTEGER            *LDB,
                    DOUBLE_COMPLEX           *X,
                    const INTEGER            *LDX,
                    DOUBLE                   *FERR,
                    DOUBLE                   *BERR,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
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
    } else if (*LDAF<std::max(INTEGER(1), *N)) {
        *INFO = -7;
    } else if (*LDB<std::max(INTEGER(1), *N)) {
        *INFO = -10;
    } else if (*LDX<std::max(INTEGER(1), *N)) {
        *INFO = -12;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZGERFS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto zA   = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    const auto zAF  = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(AF);
    const auto zB   = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(B);
    auto zX         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(X);
    auto zWORK      = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    Transpose              trans  = convertTo<Transpose>(*TRANS);
    ZConstGeMatrixView     A_     = ZConstFSView(*N, *N, zA, *LDA);
    ZConstGeMatrixView     AF_    = ZConstFSView(*N, *N, zAF, *LDAF);
    IConstDenseVectorView  IPIV_  = IConstArrayView(*N, IPIV, 1);
    ZConstGeMatrixView     B_     = ZConstFSView(*N, *NRHS, zB, *LDB);
    ZGeMatrixView          X_     = ZFSView(*N, *NRHS, zX, *LDX);
    DDenseVectorView       FERR_  = DArrayView(*NRHS, FERR, 1);
    DDenseVectorView       BERR_  = DArrayView(*NRHS, BERR, 1);
    ZDenseVectorView       WORK_  = ZArrayView(*N*2, zWORK, 1);
    DDenseVectorView       RWORK_ = DArrayView(*N, RWORK, 1);

    rfs(trans, A_, AF_, IPIV_, B_, X_, FERR_, BERR_, WORK_, RWORK_);
}


} // extern "C"

} } // namespace lapack, flens
