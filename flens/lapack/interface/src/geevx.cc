#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeev ---------------------------------------------------------------------
void
LAPACK_DECL(dgeevx)(const char       *BALANC,
                    const char       *JOBVL,
                    const char       *JOBVR,
                    const char       *_SENSE,
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
    const bool wantSNN = (*_SENSE=='N');
    const bool wantSNE = (*_SENSE=='E');
    const bool wantSNV = (*_SENSE=='V');
    const bool wantSNB = (*_SENSE=='B');

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
    SENSE::Sense      sense   = SENSE::Sense(*_SENSE);

    ASSERT(char(balance)==*BALANC);
    ASSERT(char(sense)==*_SENSE);

    typedef typename DGeMatrixView::IndexType   IndexType;

    DGeMatrixView       _A      = DFSView(*N, *N, A, *LDA);
    DDenseVectorView    _WR     = DArrayView(*N, WR, 1);
    DDenseVectorView    _WI     = DArrayView(*N, WI, 1);
    DGeMatrixView       _VL     = DFSView(*N, *N, VL, *LDVL);
    DGeMatrixView       _VR     = DFSView(*N, *N, VR, *LDVR);
    IndexType           _ILO    = *ILO;
    IndexType           _IHI    = *IHI;
    DDenseVectorView    _SCALE  = DArrayView(*N, SCALE, 1);
    DDenseVectorView    _RCONDE = DArrayView(*N, RCONDE, 1);
    DDenseVectorView    _RCONDV = DArrayView(*N, RCONDV, 1);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    IDenseVector        _IWORK(*N * 2 - 2);
    for (int i=1; i<=_IWORK.length(); ++i) {
        _IWORK(i) = IWORK[i-1];
    }

//
//  Test if work has at least minimal worksize
//
    auto ws = evx_wsq(wantVL, wantVR, sense, _A);

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

    evx(balance, wantVL, wantVR, sense, _A, _WR, _WI, _VL, _VR, _ILO, _IHI,
        _SCALE, *ABNRM, _RCONDE, _RCONDV, _WORK, _IWORK);

    for (int i=1; i<=_IWORK.length(); ++i) {
        IWORK[i-1] = _IWORK(i);
    }

    *ILO = _ILO;
    *IHI = _IHI;
}

} // extern "C"

} } // namespace lapack, flens
