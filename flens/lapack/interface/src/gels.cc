//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dgels ---------------------------------------------------------------------
void
LAPACK_DECL(dgels)(const char           *TRANS,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *B,
                   const INTEGER        *LDB,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("dgels");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lQuery = (*LWORK==0);
    const INTEGER mn = std::min(*M, *N);

    *INFO = 0;
    if (*TRANS!='N' && *TRANS!='T' && *TRANS!='C') {
        *INFO = -1;
    } else if (*M<0) {
        *INFO = -2;
    } else if (*N<0) {
        *INFO = -3;
    } else if (*NRHS<0) {
        *INFO = -4;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -6;
    } else if (*LDB<flens::max(INTEGER(1), *M, *N)) {
        *INFO = -8;
    } else if (*LWORK<std::max(INTEGER(1), mn+std::max(mn, *NRHS)) && !lQuery) {
        *INFO = -10;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGELS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    Transpose           trans  = convertTo<Transpose>(*TRANS);
    DGeMatrixView       _A     = DFSView(*M, *N, A, *LDA);

    const INTEGER numRowsB = std::max(*M,*N);
    DGeMatrixView       _B     = DFSView(numRowsB, *NRHS, B, *LDB);
    DDenseVectorView    _WORK  = DArrayView(*LWORK, WORK, 1);

    ls(trans, _A, _B, _WORK);

}

//-- zgels ---------------------------------------------------------------------
void
LAPACK_DECL(zgels)(const char           *TRANS,
                   const INTEGER        *M,
                   const INTEGER        *N,
                   const INTEGER        *NRHS,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *B,
                   const INTEGER        *LDB,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("zgels");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lQuery = (*LWORK==0);
    const INTEGER mn = std::min(*M, *N);

    *INFO = 0;
    if (*TRANS!='N' && *TRANS!='T' && *TRANS!='C') {
        *INFO = -1;
    } else if (*M<0) {
        *INFO = -2;
    } else if (*N<0) {
        *INFO = -3;
    } else if (*NRHS<0) {
        *INFO = -4;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -6;
    } else if (*LDB<flens::max(INTEGER(1), *M, *N)) {
        *INFO = -8;
    } else if (*LWORK<std::max(INTEGER(1), mn+std::max(mn, *NRHS)) && !lQuery) {
        *INFO = -10;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZGELS", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    auto zA         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zB         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(B);
    auto zWORK      = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    Transpose           trans  = convertTo<Transpose>(*TRANS);
    ZGeMatrixView       _A     = ZFSView(*M, *N, zA, *LDA);

    const INTEGER numRowsB = std::max(*M,*N);
    ZGeMatrixView       _B     = ZFSView(numRowsB, *NRHS, zB, *LDB);
    ZDenseVectorView    _WORK  = ZArrayView(*LWORK, zWORK, 1);

    ls(trans, _A, _B, _WORK);

}


} // extern "C"

} } // namespace lapack, flens
