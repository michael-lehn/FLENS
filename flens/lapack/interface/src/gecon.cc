#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgecon --------------------------------------------------------------------
void
LAPACK_DECL(dgecon)(const char       *NORM,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dgecon");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*NORM!='1' && *NORM!='O' && *NORM!='I') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    } else if (*ANORM<0) {
        *INFO = -5;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGECON", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    Norm                norm  = Norm(*NORM!='1' ? *NORM : 'O');
    DConstGeMatrixView  _A    = DConstFSView(*N, *N, *LDA, A);
    DDenseVectorView    work  = DArrayView(*N*4, WORK, 1);
    IDenseVectorView    iwork = IArrayView(*N, IWORK, 1);

    con(norm, _A, *ANORM, *RCOND, work, iwork);
}

//-- zgecon --------------------------------------------------------------------
void
LAPACK_DECL(zgecon)(const char               *NORM,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: zgecon");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*NORM!='1' && *NORM!='O' && *NORM!='I') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    } else if (*ANORM<0) {
        *INFO = -5;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZGECON", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    Norm                norm  = Norm(*NORM!='1' ? *NORM : 'O');

    const auto *zA = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto *zWORK = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZConstGeMatrixView  _A    = ZConstFSView(*N, *N, *LDA, zA);
    ZDenseVectorView    work  = ZArrayView(*N*2, zWORK, 1);
    DDenseVectorView    rwork = DArrayView(*N*2, RWORK, 1);

    con(norm, _A, *ANORM, *RCOND, work, rwork);
}

} // extern "C"

} } // namespace lapack, flens
