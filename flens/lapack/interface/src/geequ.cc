#define STR(x)      #x
#define STRING(x)   STR(x)

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

    DConstGeMatrixView  A_  = DConstFSView(*M, *N, A, *LDA);
    DDenseVectorView    R_  = DArrayView(*M, R, 1);
    DDenseVectorView    C_  = DArrayView(*N, C, 1);

    *INFO = equ(A_, R_, C_, *ROWCND, *COLCND, *AMAX);
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
    ZConstGeMatrixView  A_  = ZConstFSView(*M, *N, zA, *LDA);

    DDenseVectorView    R_  = DArrayView(*M, R, 1);
    DDenseVectorView    C_  = DArrayView(*N, C, 1);

    *INFO = equ(A_, R_, C_, *ROWCND, *COLCND, *AMAX);
}
*/

} // extern "C"

} } // namespace lapack, flens
