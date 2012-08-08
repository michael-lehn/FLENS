#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dgelsy --------------------------------------------------------------------
void
LAPACK_DECL(dgelsy)(const INTEGER    *M,
                    const INTEGER    *N,
                    const INTEGER    *NRHS,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *B,
                    const INTEGER    *LDB,
                    INTEGER          *JPVT,
                    const DOUBLE     *RCOND,
                    INTEGER          *RANK,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    std::cerr << "dgelsy" << std::endl;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    const bool lQuery = (*LWORK==0);
    const INTEGER mn = std::min(*M, *N);

    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*NRHS<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -5;
    } else if (*LDB<flens::max(INTEGER(1), *M, *N)) {
        *INFO = -7;
    }
    if (*INFO==0) {
        INTEGER lWorkMin;
        if (mn==0 || *NRHS==0) {
            lWorkMin = 1;
        } else {
            lWorkMin = mn + std::max(2*mn,
                                     std::max(*N+1, mn + *NRHS));
        }
        if (*LWORK<lWorkMin && !lQuery) {
            *INFO = -12;
        }
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGELSY", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    DGeMatrixView       _A     = DFSView(*M, *N, A, *LDA);
    const INTEGER numRowsB = std::max(*M,*N);
    DGeMatrixView       _B     = DFSView(numRowsB, *NRHS, B, *LDB);
    IDenseVectorView    _JPVT  = IArrayView(*N, JPVT, 1);
    DDenseVectorView    _WORK  = DArrayView(*LWORK, WORK, 1);

    *RANK = lsy(_A, _B, _JPVT, *RCOND, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
