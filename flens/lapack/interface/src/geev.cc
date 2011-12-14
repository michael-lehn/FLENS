#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeev ---------------------------------------------------------------------
void
LAPACK_DECL(dgeev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE               *A,
                   const INTEGER        *LDA,
                   DOUBLE               *WR,
                   DOUBLE               *WI,
                   DOUBLE               *VL,
                   const INTEGER        *LDVL,
                   DOUBLE               *VR,
                   const INTEGER        *LDVR,
                   DOUBLE               *WORK,
                   const INTEGER        *LWORK,
                   INTEGER              *INFO)
{
    DEBUG_FLENS_LAPACK("dgeev");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    bool lQuery = (*LWORK==-1);
    bool wantVL = (*JOBVL=='V');
    bool wantVR = (*JOBVR=='V');

    if ((!wantVL) && (*JOBVL!='N')) {
        *INFO = 1;
    } else if ((!wantVR) && (*JOBVR!='N')) {
        *INFO = 2;
    } else if (*N<0) {
        *INFO = 3;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 5;
    } else if (*LDVL<1 || (wantVL && *LDVL<*N)) {
        *INFO = 9;
    } else if (*LDVR<1 || (wantVR && *LDVR<*N)) {
        *INFO = 11;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEEV", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    DGeMatrixView       _A      = DFSView(*N, *N, A, *LDA);
    DDenseVectorView    _WR     = DArrayView(*N, WR, 1);
    DDenseVectorView    _WI     = DArrayView(*N, WI, 1);
    DGeMatrixView       _VL     = DFSView(*N, *N, VL, *LDVL);
    DGeMatrixView       _VR     = DFSView(*N, *N, VR, *LDVR);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

//
//  Test if work has at least minimal worksize
//
    auto ws = ev_wsq(wantVL, wantVR, _A);

    if (*LWORK<ws.first && !lQuery) {
        *INFO = 13;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEEV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    ev(wantVL, wantVR, _A, _WR, _WI, _VL, _VR, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
