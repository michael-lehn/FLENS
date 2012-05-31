#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeequ --------------------------------------------------------------------
void
LAPACK_DECL(dgeequ)(const INTEGER    *M,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *ROWCND,
                    DOUBLE           *COLCND,
                    DOUBLE           *AMAX,
                    INTEGER          *INFO)
{
//
//  Call FLENS implementation
//
    *INFO = 0;
    if (*M<0) {
        *INFO = 1;
    } else if (*N<0) {
        *INFO = 2;
    } else if (*LDA<std::max(INTEGER(1),*M)) {
        *INFO = 4;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEEQU", INFO);
        *INFO = -(*INFO);
        return;
    }

    DConstGeMatrixView  _A  = DConstFSView(*M, *N, A, *LDA);
    DDenseVectorView    _R  = DArrayView(*M, R, 1);
    DDenseVectorView    _C  = DArrayView(*N, C, 1);

    *INFO = equ(_A, _R, _C, *ROWCND, *COLCND, *AMAX);
}

//-- zgeequ --------------------------------------------------------------------
/*
void
LAPACK_DECL(zgeequ)(const INTEGER            *M,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    DOUBLE                   *R,
                    DOUBLE                   *C,
                    DOUBLE                   *ROWCND,
                    DOUBLE                   *COLCND,
                    DOUBLE                   *AMAX,
                    INTEGER                  *INFO)
{
//
//  Call FLENS implementation
//
    *INFO = 0;
    if (*M<0) {
        *INFO = 1;
    } else if (*N<0) {
        *INFO = 2;
    } else if (*LDA<std::max(INTEGER(1),*M)) {
        *INFO = 4;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEEQU", INFO);
        *INFO = -(*INFO);
        return;
    }

    const auto *zA = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    ZConstGeMatrixView  _A  = ZConstFSView(*M, *N, zA, *LDA);

    DDenseVectorView    _R  = DArrayView(*M, R, 1);
    DDenseVectorView    _C  = DArrayView(*N, C, 1);

    *INFO = equ(_A, _R, _C, *ROWCND, *COLCND, *AMAX);
}
*/

} // extern "C"

} } // namespace lapack, flens
