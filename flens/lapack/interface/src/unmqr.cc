#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- zunmqr --------------------------------------------------------------------
void
LAPACK_DECL(zunmqr)(const char               *SIDE,
                    const char               *TRANS,
                    const INTEGER            *M,
                    const INTEGER            *N,
                    const INTEGER            *K,
                    DOUBLE_COMPLEX           *A,
                    const INTEGER            *LDA,
                    const DOUBLE_COMPLEX     *TAU,
                    DOUBLE_COMPLEX           *C,
                    const INTEGER            *LDC,
                    DOUBLE_COMPLEX           *WORK,
                    const INTEGER            *LWORK,
                    INTEGER                  *INFO)
{
    DEBUG_FLENS_LAPACK("zunmqr");

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
    } else if (!noTran && *TRANS!='C') {
        *INFO = -2;
    } else if (*M<0) {
        *INFO = -3;
    } else if (*N<0) {
        *INFO = -4;
    } else if (*K<0 || *K>nq) {
        *INFO = -5;
    } else if (*LDA<max(INTEGER(1), nq)) {
        *INFO = -7;
    } else if (*LDC<max(INTEGER(1), *M)) {
        *INFO = -10;
    } else if ((*LWORK<max(INTEGER(1), nw)) && (!lQuery)) {
        *INFO = -12;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZUNMQR", INFO);
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
    Side      side  = (*SIDE=='L') ? Left : Right;
    Transpose trans = (*TRANS=='N') ? NoTrans : ConjTrans;
    const INTEGER numRowsA = (side==Left) ? *M : *N;

    auto zA         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    const auto zTAU = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(TAU);
    auto zC         = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(C);
    auto zWORK      = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZGeMatrixView          _A      = ZFSView(numRowsA, *K, zA, *LDA);
    ZConstDenseVectorView  _TAU    = ZConstArrayView(*K, zTAU, 1);
    ZGeMatrixView          _C      = ZFSView(*M, *N, zC, *LDC);
    ZDenseVectorView       _WORK   = ZArrayView(*LWORK, zWORK, 1);

    unmqr(side, trans, _A, _TAU, _C, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
