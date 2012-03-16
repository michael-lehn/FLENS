#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgetrf --------------------------------------------------------------------
void
LAPACK_DECL(dgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgetrf");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -4;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGETRF", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//

    DGeMatrixView       _A      = DFSView(*M, *N, A, *LDA);
    IDenseVectorView    _IPIV   = IArrayView(*M, IPIV, 1);

    *INFO = trf(_A, _IPIV);
}

//-- zgetrf --------------------------------------------------------------------
void
LAPACK_DECL(zgetrf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *IPIV,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("zgetrf");
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -4;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZGETRF", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//

    ZGeMatrixView       _A      = ZFSView(*M, *N, A, *LDA);
    IDenseVectorView    _IPIV   = IArrayView(*M, IPIV, 1);

    *INFO = trf(_A, _IPIV);
}


} // extern "C"

} } // namespace lapack, flens
