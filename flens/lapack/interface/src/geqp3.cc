#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgeqp3");

    using std::max;
    using std::min;
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
    } else if (*LWORK<3*(*N)+1) {
        *INFO = -8;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGEQP3", INFO);
        *INFO = -(*INFO);
        return;
    }

    DGeMatrixView       _A      = DFSView(*M, *N, A, *LDA);
    IDenseVectorView    _JPVT   = IArrayView(*N, JPVT, 1);
    DDenseVectorView    _TAU    = DArrayView(min(*M,*N), TAU, 1);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    qp3(_A, _JPVT, _TAU, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
