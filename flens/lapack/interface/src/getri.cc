#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgetri --------------------------------------------------------------------
void
LAPACK_DECL(dgetri)(const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const INTEGER    *IPIV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgetri");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    bool lQuery = (*LWORK==-1);

    if (*N<0) {
        *INFO = 1;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 3;
    } else if (*LWORK<max(INTEGER(1),*N) && !lQuery) {
        *INFO = 6;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGETRI", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    DGeMatrixView          _A      = DFSView(*N, *N, A, *LDA);
    IConstDenseVectorView  _IPIV   = IConstArrayView(*N, IPIV, 1);
    DDenseVectorView       _WORK   = DArrayView(*LWORK, WORK, 1);

//
//  Call FLENS implementation
//
    tri(_A, _IPIV, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
