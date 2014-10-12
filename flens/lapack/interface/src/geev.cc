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
    DGeMatrixView       A_      = DFSView(*N, *N, A, *LDA);
    DDenseVectorView    WR_     = DArrayView(*N, WR, 1);
    DDenseVectorView    WI_     = DArrayView(*N, WI, 1);
    DGeMatrixView       VL_     = DFSView(*N, *N, VL, *LDVL);
    DGeMatrixView       VR_     = DFSView(*N, *N, VR, *LDVR);
    DDenseVectorView    WORK_   = DArrayView(*LWORK, WORK, 1);

//
//  Test if work has at least minimal worksize
//
    auto ws = ev_wsq(wantVL, wantVR, A_);

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
    ev(wantVL, wantVR, A_, WR_, WI_, VL_, VR_, WORK_);
}

//-- zgeev ---------------------------------------------------------------------
void
LAPACK_DECL(zgeev)(const char           *JOBVL,
                   const char           *JOBVR,
                   const INTEGER        *N,
                   DOUBLE_COMPLEX       *A,
                   const INTEGER        *LDA,
                   DOUBLE_COMPLEX       *W,
                   DOUBLE_COMPLEX       *VL,
                   const INTEGER        *LDVL,
                   DOUBLE_COMPLEX       *VR,
                   const INTEGER        *LDVR,
                   DOUBLE_COMPLEX       *WORK,
                   const INTEGER        *LWORK,
                   DOUBLE               *RWORK,
                   INTEGER              *INFO)
{
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
        *INFO = 8;
    } else if (*LDVR<1 || (wantVR && *LDVR<*N)) {
        *INFO = 10;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEEV", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zW     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(W);
    auto zVL    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(VL);
    auto zVR    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(VR);
    auto zWORK  = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZGeMatrixView       A_      = ZFSView(*N, *N, zA, *LDA);
    ZDenseVectorView    W_      = ZArrayView(*N, zW, 1);
    ZGeMatrixView       VL_     = ZFSView(*N, *N, zVL, *LDVL);
    ZGeMatrixView       VR_     = ZFSView(*N, *N, zVR, *LDVR);
    ZDenseVectorView    WORK_   = ZArrayView(*LWORK, zWORK, 1);
    DDenseVectorView    RWORK_  = DArrayView(2*(*N), RWORK, 1);

//
//  Test if work has at least minimal worksize
//
    auto ws = ev_wsq(wantVL, wantVR, A_);

    if (*LWORK<ws.first && !lQuery) {
        *INFO = 12;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEEV", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    ev(wantVL, wantVR, A_, W_, VL_, VR_, WORK_, RWORK_);
}

} // extern "C"

} } // namespace lapack, flens
