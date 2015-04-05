#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dorglq --------------------------------------------------------------------
void
LAPACK_DECL(dorglq)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dorglq");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lQuery = (*LWORK==-1);

    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<*M) {
        *INFO = -2;
    } else if (*K<0 || *K>*M) {
        *INFO = -3;
    } else if (*LDA<max(INTEGER(1), *M)) {
        *INFO = -5;
    } else if (*LWORK<max(INTEGER(1), *M) && !lQuery) {
        *INFO = -8;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DORGLQ", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Handle worksize query
//
    if (lQuery) {
        // TODO: implement wsq
        ASSERT(0);
    }
//
//  Call FLENS implementation
//
    DGeMatrixView          _A      = DFSView(*M, *N, *LDA, A);
    DConstDenseVectorView  _TAU    = DConstArrayView(*K, TAU, 1);
    DDenseVectorView       _WORK   = DArrayView(*LWORK, WORK, 1);

    orglq(_A, _TAU, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
