#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dormlq --------------------------------------------------------------------
void
LAPACK_DECL(dormlq)(const char       *SIDE,
                    const char       *TRANS,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *K,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *TAU,
                    DOUBLE           *C,
                    const INTEGER    *LDC,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dormlq");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool left   = (*SIDE=='L');
    const bool noTran = (*TRANS=='N');
    const bool lQuery = (*LWORK==-1);

    INTEGER nq, nw;

    if (left) {
        nq = *M;
        nw = *N;
    } else {
        nq = *N;
        nw = *M;
    }
    if (!left && *SIDE!='R') {
        *INFO = -1;
    } else if (!noTran && *TRANS!='T') {
        *INFO = -2;
    } else if (*M<0) {
        *INFO = -3;
    } else if (*N<0) {
        *INFO = -4;
    } else if (*K<0 || *K>nq) {
        *INFO = -5;
    } else if (*LDA<max(INTEGER(1), *K)) {
        *INFO = -7;
    } else if (*LDC<max(INTEGER(1), *M)) {
        *INFO = -10;
    } else if ((*LWORK<max(INTEGER(1), nw)) && (!lQuery)) {
        *INFO = -12;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DORMLQ", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Handle worksize query
//
    if (lQuery) {
        // TODO: implement lqf_wsq
        ASSERT(0);
    }
//
//  Call FLENS implementation
//
    Side      side  = (*SIDE=='L') ? Left : Right;
    Transpose trans = (*TRANS=='N') ? NoTrans : Trans;
    const INTEGER numColsA = (side==Left) ? *M : *N;

    DGeMatrixView          _A      = DFSView(*K, numColsA, A, *LDA);
    DConstDenseVectorView  _TAU    = DConstArrayView(*K, TAU, 1);
    DGeMatrixView          _C      = DFSView(*M, *N, C, *LDC);
    DDenseVectorView       _WORK   = DArrayView(*LWORK, WORK, 1);

    ormlq(side, trans, _A, _TAU, _C, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
