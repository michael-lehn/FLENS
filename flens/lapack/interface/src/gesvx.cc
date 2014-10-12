#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgesvx --------------------------------------------------------------------
void
LAPACK_DECL(dgesvx)(const char       *FACT,
                    const char       *TRANS,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *AF,
                    const INTEGER    *LDAF,
                    INTEGER          *IPIV,
                    char             *EQUED,
                    DOUBLE           *R,
                    DOUBLE           *C,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    DOUBLE           *X,
                    const INTEGER    *LDX,
                    DOUBLE           *RCOND,
                    DOUBLE           *FERR,
                    DOUBLE           *BERR,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
//
//  Test the input parameters so that we pass LAPACK error checks
//
    typedef DOUBLE   T;
    typedef INTEGER  IndexType;

    const T  Zero(0), One(1);
    const T bigNum = One / lamch<T>(SafeMin);

    bool rowEqu, colEqu;

    if (*FACT=='N' || *FACT=='E') {
        *EQUED = 'N';
        rowEqu = false;
        colEqu = false;
    } else {
        rowEqu = (*EQUED=='R' || *EQUED=='B');
        colEqu = (*EQUED=='C' || *EQUED=='B');
    }

    *INFO = 0;
    if (*FACT!='F' && *FACT!='N' && *FACT!='E') {
        *INFO = -1;
    } else if (*TRANS!='N' && *TRANS!='T' && *TRANS!='C') {
        *INFO = -2;
    } else if (*N<0) {
        *INFO = -3;
    } else if (*NRHS<0) {
        *INFO = -4;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -6;
    } else if (*LDAF<std::max(INTEGER(1), *N)) {
        *INFO = -8;
    } else if (*FACT=='F' && !(rowEqu || colEqu || *EQUED=='N')) {
        *INFO = -10;
    } else {
        if (rowEqu) {
            T rcMin = bigNum;
            for (IndexType j=0; j<*N; ++j) {
                rcMin = min(rcMin, R[j]);
            }
            if (rcMin<=Zero) {
                *INFO = -11;
            }
        }
        if (colEqu && *INFO==0) {
            T rcMin = bigNum;
            for (IndexType j=0; j<*N; ++j) {
                rcMin = min(rcMin, C[j]);
            }
            if (rcMin<=Zero) {
                *INFO = -12;
            }
        }
        if (*INFO==0) {
            if (*LDB<std::max(INTEGER(1), *N)) {
                *INFO = -14;
            } else if (*LDX<std::max(INTEGER(1), *N)) {
                *INFO = -16;
            }
        }
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGESVX", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    SVX::Fact          FACT_  = SVX::Fact(*FACT);
    Transpose          TRANS_ = cxxblas::getCxxBlasEnum<Transpose>(*TRANS);
    DGeMatrixView      A_     = DFSView(*N, *N, A, *LDA);
    DGeMatrixView      AF_    = DFSView(*N, *N, AF, *LDAF);
    IDenseVectorView   IPIV_  = IArrayView(*N, IPIV, 1);
    SVX::Equilibration EQUED_ = SVX::Equilibration(*EQUED);
    DDenseVectorView   R_     = DArrayView(*N, R, 1);
    DDenseVectorView   C_     = DArrayView(*N, C, 1);
    DGeMatrixView      B_     = DFSView(*N, *NRHS, B, *LDB);
    DGeMatrixView      X_     = DFSView(*N, *NRHS, X, *LDX);
    DDenseVectorView   FERR_  = DArrayView(*NRHS, FERR, 1);
    DDenseVectorView   BERR_  = DArrayView(*NRHS, BERR, 1);
    DDenseVectorView   WORK_  = DArrayView(*N*4, WORK, 1);
    IDenseVectorView   IWORK_ = IArrayView(*N, IWORK, 1);

    *INFO = svx(FACT_, TRANS_, A_, AF_, IPIV_, EQUED_, R_, C_, B_, X_, *RCOND,
                FERR_, BERR_, WORK_, IWORK_);

    *EQUED = char(EQUED_);
}

} // extern "C"

} } // namespace lapack, flens
