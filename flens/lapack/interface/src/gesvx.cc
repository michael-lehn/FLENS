//#define CXXBLAS_DEBUG_OUT(x)      std::cerr << x << std::endl;

#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

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
    DEBUG_FLENS_LAPACK("dgesvx");
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
    SVX::Fact          _FACT  = SVX::Fact(*FACT);
    Transpose          _TRANS = cxxblas::getCxxBlasEnum<Transpose>(*TRANS);
    DGeMatrixView      _A     = DFSView(*N, *N, A, *LDA);
    DGeMatrixView      _AF    = DFSView(*N, *N, AF, *LDAF);
    IDenseVectorView   _IPIV  = IArrayView(*N, IPIV, 1);
    SVX::Equilibration _EQUED = SVX::Equilibration(*EQUED);
    DDenseVectorView   _R     = DArrayView(*N, R, 1);
    DDenseVectorView   _C     = DArrayView(*N, C, 1);
    DGeMatrixView      _B     = DFSView(*N, *NRHS, B, *LDB);
    DGeMatrixView      _X     = DFSView(*N, *NRHS, X, *LDX);
    DDenseVectorView   _FERR  = DArrayView(*NRHS, FERR, 1);
    DDenseVectorView   _BERR  = DArrayView(*NRHS, BERR, 1);
    DDenseVectorView   _WORK  = DArrayView(*N*4, WORK, 1);
    IDenseVectorView   _IWORK = IArrayView(*N, IWORK, 1);

    *INFO = svx(_FACT, _TRANS, _A, _AF, _IPIV, _EQUED, _R, _C, _B, _X, *RCOND,
                _FERR, _BERR, _WORK, _IWORK);

    *EQUED = char(_EQUED);
}

} // extern "C"

} } // namespace lapack, flens
