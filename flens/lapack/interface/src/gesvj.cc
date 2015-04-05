#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgetrs --------------------------------------------------------------------
void
LAPACK_DECL(dgesvj)(const char       *JOBA,
                    const char       *JOBU,
                    const char       *JOBV,
                    const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *SVA,
                    const INTEGER    *MV,
                    DOUBLE           *V,
                    const INTEGER    *LDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dgesvj");

//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lsvec = (*JOBU=='U');
    const bool uctol = (*JOBU=='C');
    const bool rsvec = (*JOBV=='V');
    const bool applv = (*JOBV=='A');
    const bool upper = (*JOBA=='U');
    const bool lower = (*JOBA=='L');

    *INFO = 0;
    if (!(upper || lower || *JOBA=='G')) {
        *INFO = -1;
    } else if (!(lsvec || uctol || *JOBU=='N')) {
        *INFO = -2;
    } else if (!(rsvec || applv || *JOBV=='N')) {
        *INFO = -3;
    } else if (*M<0) {
        *INFO = -4;
    } else if ((*N<0) || (*N>*M)) {
        *INFO = -5;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -7;
    } else if (*MV<INTEGER(0)) {
        *INFO = -9;
    } else if ((rsvec && *LDV<*N) || (applv && *LDV<*MV)) {
        *INFO = -11;
    } else if (uctol && (*WORK<=DOUBLE(1))) {
        *INFO = -12;
    } else if (*LWORK<std::max(*M+*N,INTEGER(6))) {
        *INFO = -13;
    }

    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGESVJ", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Call FLENS implementation
//
    INTEGER mv = 0;

    if (*JOBV=='A') {
        mv = *MV;
    }
    if (*JOBV=='V') {
        mv = *N;
    }

    SVJ::TypeA        typeA  = SVJ::TypeA(*JOBA);
    SVJ::JobU         jobU   = SVJ::JobU(*JOBU);
    SVJ::JobV         jobV   = SVJ::JobV(*JOBV);
    DGeMatrixView     _A     = DFSView(*M, *N, *LDA, A);
    DDenseVectorView  _sva   = DArrayView(*N, SVA, INTEGER(1));
    DGeMatrixView     _V     = DFSView(mv, *N, *LDV, V);
    DDenseVectorView  _work  = DArrayView(*LWORK, WORK, INTEGER(1));

    *INFO = svj(typeA, jobU, jobV, _A, _sva, _V, _work);

    if (*INFO<0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGESVJ", INFO);
        *INFO = -(*INFO);
        return;
    }
}

} // extern "C"

} } // namespace lapack, flens
