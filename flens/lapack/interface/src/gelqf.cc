#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgelqf --------------------------------------------------------------------
void
LAPACK_DECL(dgelqf)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    DEBUG_FLENS_LAPACK("dgelqf");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    bool lQuery = (*LWORK==-1);

    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -4;
    } else if ((*LWORK<max(INTEGER(1), *M)) && (!lQuery)) {
        *INFO = -7;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGELQF", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Handle worksize query
//
    if (lQuery) {
        // TODO: implement lqf_wsq
        ASSERT(0);
    }
//
//  Call FLENS implementation
//
    DGeMatrixView       _A      = DFSView(*M, *N, A, *LDA);
    DDenseVectorView    _TAU    = DArrayView(min(*M,*N), TAU, 1);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    lqf(_A, _TAU, _WORK);
}

} // extern "C"

} } // namespace lapack, flens
