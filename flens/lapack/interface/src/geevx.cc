#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeevx --------------------------------------------------------------------
void
LAPACK_DECL(dgeevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE_,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    DOUBLE           *VL,
                    const INTEGER    *LDVL,
                    DOUBLE           *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool lQuery = (*LWORK==-1);
    const bool wantVL = (*JOBVL=='V');
    const bool wantVR = (*JOBVR=='V');
    const bool wantSNN = (*SENSE_=='N');
    const bool wantSNE = (*SENSE_=='E');
    const bool wantSNV = (*SENSE_=='V');
    const bool wantSNB = (*SENSE_=='B');

    if (*BALANC!='N' && *BALANC!='S' && *BALANC!='P' && *BALANC!='B') {
        *INFO = 1;
    } else if ((!wantVL) && (*JOBVL!='N')) {
        *INFO = 2;
    } else if ((!wantVR) && (*JOBVR!='N')) {
        *INFO = 3;
    } else if (!(wantSNN || wantSNE || wantSNB || wantSNV)
            || ((wantSNE || wantSNB ) && !(wantVL && wantVR)))
    {
        *INFO = 4;
    } else if (*N<0) {
        *INFO = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 7;
    } else if (*LDVL<1 || (wantVL && *LDVL<*N)) {
        *INFO = 11;
    } else if (*LDVR<1 || (wantVR && *LDVR<*N)) {
        *INFO = 13;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEEVX", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    BALANCE::Balance  balance = BALANCE::Balance(*BALANC);
    SENSE::Sense      sense   = SENSE::Sense(*SENSE_);

    ASSERT(char(balance)==*BALANC);
    ASSERT(char(sense)==*SENSE_);

    typedef typename DGeMatrixView::IndexType   IndexType;

    DGeMatrixView       A_      = DFSView(*N, *N, A, *LDA);
    DDenseVectorView    WR_     = DArrayView(*N, WR, 1);
    DDenseVectorView    WI_     = DArrayView(*N, WI, 1);
    DGeMatrixView       VL_     = DFSView(*N, *N, VL, *LDVL);
    DGeMatrixView       VR_     = DFSView(*N, *N, VR, *LDVR);
    IndexType           ILO_    = *ILO;
    IndexType           IHI_    = *IHI;
    DDenseVectorView    SCALE_  = DArrayView(*N, SCALE, 1);
    DDenseVectorView    RCONDE_ = DArrayView(*N, RCONDE, 1);
    DDenseVectorView    RCONDV_ = DArrayView(*N, RCONDV, 1);
    DDenseVectorView    WORK_   = DArrayView(*LWORK, WORK, 1);

    IDenseVector        IWORK_(*N * 2 - 2);
    for (int i=1; i<=IWORK_.length(); ++i) {
        IWORK_(i) = IWORK[i-1];
    }

//
//  Test if work has at least minimal worksize
//
    auto ws = evx_wsq(wantVL, wantVR, sense, A_);

    if (*LWORK<ws.first && !lQuery) {
        *INFO = 21;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEEVX", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//

    evx(balance, wantVL, wantVR, sense, A_, WR_, WI_, VL_, VR_, ILO_, IHI_,
        SCALE_, *ABNRM, RCONDE_, RCONDV_, WORK_, IWORK_);

    for (int i=1; i<=IWORK_.length(); ++i) {
        IWORK[i-1] = IWORK_(i);
    }

    *ILO = ILO_;
    *IHI = IHI_;
}

//-- zgeevx --------------------------------------------------------------------
void
LAPACK_DECL(zgeevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *SENSE_,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    DOUBLE_COMPLEX   *W,
                    DOUBLE_COMPLEX   *VL,
                    const INTEGER    *LDVL,
                    DOUBLE_COMPLEX   *VR,
                    const INTEGER    *LDVR,
                    INTEGER          *ILO,
                    INTEGER          *IHI,
                    DOUBLE           *SCALE,
                    DOUBLE           *ABNRM,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO)
{
    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool lQuery = (*LWORK==-1);
    const bool wantVL = (*JOBVL=='V');
    const bool wantVR = (*JOBVR=='V');
    const bool wantSNN = (*SENSE_=='N');
    const bool wantSNE = (*SENSE_=='E');
    const bool wantSNV = (*SENSE_=='V');
    const bool wantSNB = (*SENSE_=='B');

    if (*BALANC!='N' && *BALANC!='S' && *BALANC!='P' && *BALANC!='B') {
        *INFO = 1;
    } else if ((!wantVL) && (*JOBVL!='N')) {
        *INFO = 2;
    } else if ((!wantVR) && (*JOBVR!='N')) {
        *INFO = 3;
    } else if (!(wantSNN || wantSNE || wantSNB || wantSNV)
            || ((wantSNE || wantSNB ) && !(wantVL && wantVR)))
    {
        *INFO = 4;
    } else if (*N<0) {
        *INFO = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 7;
    } else if (*LDVL<1 || (wantVL && *LDVL<*N)) {
        *INFO = 10;
    } else if (*LDVR<1 || (wantVR && *LDVR<*N)) {
        *INFO = 12;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEEVX", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    BALANCE::Balance  balance = BALANCE::Balance(*BALANC);
    SENSE::Sense      sense   = SENSE::Sense(*SENSE_);

    ASSERT(char(balance)==*BALANC);
    ASSERT(char(sense)==*SENSE_);

    typedef typename DGeMatrixView::IndexType   IndexType;

    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zW     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(W);
    auto zVL    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(VL);
    auto zVR    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(VR);
    auto zWORK  = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZGeMatrixView       A_      = ZFSView(*N, *N, zA, *LDA);
    ZDenseVectorView    W_      = ZArrayView(*N, zW, 1);
    ZGeMatrixView       VL_     = ZFSView(*N, *N, zVL, *LDVL);
    ZGeMatrixView       VR_     = ZFSView(*N, *N, zVR, *LDVR);
    IndexType           ILO_    = *ILO;
    IndexType           IHI_    = *IHI;
    DDenseVectorView    SCALE_  = DArrayView(*N, SCALE, 1);
    DDenseVectorView    RCONDE_ = DArrayView(*N, RCONDE, 1);
    DDenseVectorView    RCONDV_ = DArrayView(*N, RCONDV, 1);
    ZDenseVectorView    WORK_   = ZArrayView(*LWORK, zWORK, 1);
    DDenseVectorView    RWORK_  = DArrayView(2*(*N), RWORK, 1);

//
//  Test if work has at least minimal worksize
//
    auto ws = evx_wsq(wantVL, wantVR, sense, A_);

    if (*LWORK<ws.first && !lQuery) {
        *INFO = 20;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEEVX", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//

    evx(balance, wantVL, wantVR, sense, A_, W_, VL_, VR_, ILO_, IHI_,
        SCALE_, *ABNRM, RCONDE_, RCONDV_, WORK_, RWORK_);

    *ILO = ILO_;
    *IHI = IHI_;
}

} // extern "C"

} } // namespace lapack, flens
