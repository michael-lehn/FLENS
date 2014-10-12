#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgejsv --------------------------------------------------------------------
void
LAPACK_DECL(dgejsv)(const char       *JOBA,
                    const char       *JOBU,
                    const char       *JOBV,
                    const char       *JOBR,
                    const char       *JOBT,
                    const char       *JOBP,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *SVA,
                    DOUBLE           *U,
                    const INTEGER    *LDU,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lsvec  = (*JOBU=='U') || (*JOBU=='F');
    const bool jracc  = (*JOBV=='J');
    const bool rsvec  = (*JOBV=='V') || jracc;
    const bool rowpiv = (*JOBA=='F') || (*JOBA=='G');
    const bool l2rank = (*JOBA=='R');
    const bool l2aber = (*JOBA=='A');
    const bool errest = (*JOBA=='E') || (*JOBA=='G');
    const bool l2tran = (*JOBT=='T');
    const bool l2kill = (*JOBR=='R');
    const bool defr   = (*JOBR=='N');
    const bool l2pert = (*JOBP=='P');

    const int m = *M;
    const int n = *N;

    *INFO = 0;
    if (!(rowpiv || l2rank || l2aber || errest || *JOBA=='C')) {
        *INFO = -1;
    } else if (!(lsvec || *JOBU=='N' || *JOBU=='W')) {
        *INFO = -2;
    } else if (!(rsvec || *JOBV=='N' || *JOBV=='W') || (jracc && !lsvec)) {
        *INFO = -3;
    } else if (!(l2kill || defr)) {
        *INFO = -4;
    } else if (!(l2tran || *JOBT=='N')) {
        *INFO = -5;
    } else if (!(l2pert || *JOBP=='N')) {
        *INFO = -6;
    } else if (m<0) {
        *INFO = -7;
    } else if ((n<0) || (n>m)) {
        *INFO = -8;
    } else if (*LDA<m) {
        *INFO = -10;
    } else if (lsvec && *LDU<m) {
        *INFO = -13;
    } else if (rsvec && *LDV<n) {
        *INFO = -14;
    } else if ((!(lsvec || rsvec || errest) || (*LWORK < max(7,4*n+1,2*m+n)))
        || (!(lsvec || rsvec) || errest || (*LWORK < max(7,4*n+n*n,2*m+n)))
        || (lsvec || (!rsvec) || (*LWORK < max(7,2*m+n,4*n+1)))
        || (rsvec || (!lsvec) || (*LWORK < max(7,2*m+n,4*n+1)))
        || (lsvec || rsvec || (!jracc) || (*LWORK<max(2*m+n,6*n+2*n*n)))
        || (lsvec || rsvec || jracc || *LWORK<max(2*m+n,4*n+n*n,2*n+n*n+6)))
    {
        *INFO = -17;
    }

    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGEJSV", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Call FLENS implementation
//
    JSV::Accuracy  accuracy    = JSV::Accuracy(*JOBA);
    JSV::JobU  jobU            = JSV::JobU(*JOBU);
    JSV::JobV  jobV            = JSV::JobV(*JOBV);
    const bool restrictedRange = (*JOBR=='R');
    const bool considerTransA  = (*JOBT=='T');
    const bool perturb         = (*JOBR=='P');

    DGeMatrixView     A_       = DFSView(m, n, A, *LDA);
    DDenseVectorView  sva_     = DArrayView(n, SVA, INTEGER(1));
    DGeMatrixView     U_       = DFSView(m, n, U, *LDU);
    DGeMatrixView     V_       = DFSView(n, n, V, *LDV);
    DDenseVectorView  work_    = DArrayView(*LWORK, WORK, INTEGER(1));
    IDenseVectorView  iwork_   = IArrayView(m+3*n, IWORK, INTEGER(1));

    *INFO = jsv(accuracy, jobU, jobV,
                restrictedRange, considerTransA, perturb,
                A_, sva_, U_, V_, work_, iwork_);

    if (*INFO<0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGEJSV", INFO);
        *INFO = -(*INFO);
        return;
    }
}

} // extern "C"

} } // namespace lapack, flens
