#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(dgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE           *TAU,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dgeqp3");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -4;
    } else if (*LWORK<3*(*N)+1) {
        *INFO = -8;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DGEQP3", INFO);
        *INFO = -(*INFO);
        return;
    }

    DGeMatrixView       _A      = DFSView(*M, *N, *LDA, A);
    IDenseVectorView    _JPVT   = IArrayView(*N, JPVT, 1);
    DDenseVectorView    _TAU    = DArrayView(min(*M,*N), TAU, 1);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    qp3(_A, _JPVT, _TAU, _WORK);
}

//-- zgeqrf --------------------------------------------------------------------
void
LAPACK_DECL(zgeqp3)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *JPVT,
                    DOUBLE_COMPLEX   *TAU,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: zgeqp3");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *M)) {
        *INFO = -4;
    } else if (*LWORK<(*N)+1) {
        *INFO = -8;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZGEQP3", INFO);
        *INFO = -(*INFO);
        return;
    }

    auto zA    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zTAU  = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(TAU);
    auto zWORK = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZGeMatrixView       _A      = ZFSView(*M, *N, *LDA, zA);
    IDenseVectorView    _JPVT   = IArrayView(*N, JPVT, 1);
    ZDenseVectorView    _TAU    = ZArrayView(min(*M,*N), zTAU, 1);
    ZDenseVectorView    _WORK   = ZArrayView(*LWORK, zWORK, 1);
    DDenseVectorView    _RWORK  = DArrayView(2*(*N), RWORK, 1);

    qp3(_A, _JPVT, _TAU, _WORK, _RWORK);
}

} // extern "C"

} } // namespace lapack, flens
