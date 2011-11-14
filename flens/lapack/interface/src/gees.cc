#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

template <typename T>
struct SelectFunction
{
    typedef LOGICAL (* Function)(const T *, const T *);

    SelectFunction(Function _select)
        : select(_select)
    {
    }

    bool
    operator()(const T &a, const T &b)
    {
        return select(&a, &b);
    }

    Function  select;
};


extern "C" {

//-- dgeev ---------------------------------------------------------------------
void
LAPACK_DECL(dgees)(const char       *JOBVS,
                   const char       *SORT,
                   LOGICAL          (*SELECT)(const DOUBLE *, const DOUBLE *),
                   const INTEGER    *N,
                   DOUBLE           *A,
                   const INTEGER    *LDA,
                   INTEGER          *SDIM,
                   DOUBLE           *WR,
                   DOUBLE           *WI,
                   DOUBLE           *VS,
                   const INTEGER    *LDVS,
                   DOUBLE           *WORK,
                   const INTEGER    *LWORK,
                   LOGICAL          *BWORK,
                   INTEGER          *INFO)
{
    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    bool lQuery = (*LWORK==-1);
    bool wantVS = (*JOBVS=='V');
    bool wantST = (*SORT=='S');

    if ((!wantVS) && (*JOBVS!='N')) {
        *INFO = 1;
    } else if ((!wantST) && (*SORT!='N')) {
        *INFO = 2;
    } else if (*N<0) {
        *INFO = 4;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 6;
    } else if (*LDVS<1 || (wantVS && *LDVS<*N)) {
        *INFO = 11;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEES", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    typedef DGeMatrixView::IndexType IndexType;

    DGeMatrixView       _A      = DFSView(*N, *N, A, *LDA);
    IndexType           _SDIM   = *SDIM;
    DDenseVectorView    _WR     = DArrayView(*N, WR, 1);
    DDenseVectorView    _WI     = DArrayView(*N, WI, 1);
    DGeMatrixView       _VS     = DFSView(*N, *N, VS, *LDVS);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    BDenseVector        _BWORK(*N);
    for (INTEGER i=1; i<*N; ++i) {
        _BWORK(i) = BWORK[i-1];
    }

//
//  Test if work has at least minimal worksize
//
    auto ws = es_wsq(wantVS, _A);

    if (*LWORK<ws.first && !lQuery) {
        *INFO = 13;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEES", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    SelectFunction<DOUBLE> select(SELECT);

    es(wantVS, wantST, select, _A, _SDIM, _WR, _WI, _VS, _WORK, _BWORK);

    *SDIM = _SDIM;
    for (INTEGER i=1; i<*N; ++i) {
        BWORK[i-1] = _BWORK(i);
    }
}

} // extern "C"

} } // namespace lapack, flens
