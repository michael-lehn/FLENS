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

    LOGICAL
    operator()(const T &a, const T &b)
    {
        return select(&a, &b);
    }

    Function  select;
};

template <typename T>
struct SelectFunction<std::complex<T> >
{
    typedef std::complex<T> CT;

    typedef LOGICAL (* LapackFunction)(const T *);
    typedef LOGICAL (* Function)(const CT *);

    SelectFunction(LapackFunction _select)
        : select(reinterpret_cast<Function>(_select))
    {
    }

    LOGICAL
    operator()(const CT &a)
    {
        return select(reinterpret_cast<const T *>(&a));
    }

    Function  select;
};

extern "C" {

//-- dgeesx --------------------------------------------------------------------
void
LAPACK_DECL(dgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    LOGICAL          (*SELECT)(const DOUBLE *, const DOUBLE *),
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    INTEGER          *SDIM,
                    DOUBLE           *WR,
                    DOUBLE           *WI,
                    DOUBLE           *VS,
                    const INTEGER    *LDVS,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE           *WORK,
                    const INTEGER    *LWORK,
                    INTEGER          *IWORK,
                    const INTEGER    *LIWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dgeesx");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool lQuery = (*LWORK==-1) || (*LIWORK==-1);
    const bool wantVS = (*JOBVS=='V');
    const bool wantST = (*SORT=='S');

    const bool wantSN = (*SENSE=='N');
    const bool wantSE = (*SENSE=='E');
    const bool wantSV = (*SENSE=='V');
    const bool wantSB = (*SENSE=='B');

    if ((!wantVS) && (*JOBVS!='N')) {
        *INFO = 1;
    } else if ((!wantST) && (*SORT!='N')) {
        *INFO = 2;
    } else if (!(wantSN || wantSE || wantSV || wantSB)
            || ( !wantST && !wantSN))
    {
        *INFO = 4;
    } else if (*N<0) {
        *INFO = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 7;
    } else if (*LDVS<1 || (wantVS && *LDVS<*N)) {
        *INFO = 12;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEESX", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    typedef DGeMatrixView::IndexType IndexType;

    SENSE::Sense sense = SENSE::Sense(*SENSE);

    DGeMatrixView       _A      = DFSView(*N, *N, *LDA, A);
    IndexType           _SDIM   = *SDIM;
    DDenseVectorView    _WR     = DArrayView(*N, WR, 1);
    DDenseVectorView    _WI     = DArrayView(*N, WI, 1);
    DGeMatrixView       _VS     = DFSView(*N, *N, *LDVS, VS);
    DDenseVectorView    _WORK   = DArrayView(*LWORK, WORK, 1);

    IDenseVector        _IWORK(*LIWORK);
    for (INTEGER i=1; i<*LIWORK; ++i) {
        _IWORK(i) = IWORK[i-1];
    }

    BDenseVector        _BWORK(*N);
    for (INTEGER i=1; i<*N; ++i) {
        _BWORK(i) = BWORK[i-1];
    }

//
//  Test if work has at least minimal worksize
//
    auto wsQuery = esx_wsq(wantVS, sense, _A);
    INTEGER minWork = wsQuery.first;

    if (*LWORK<minWork && !lQuery) {
        *INFO = 16;
    } else if (*LIWORK<1 && !lQuery) {
        *INFO = 18;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("DGEESX", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    SelectFunction<DOUBLE> select(SELECT);

    esx(wantVS, wantST, select, sense, _A, _SDIM, _WR, _WI, _VS,
        *RCONDE, *RCONDV, _WORK, _IWORK, _BWORK);

    *SDIM   = _SDIM;
    for (INTEGER i=1; i<*LIWORK; ++i) {
        IWORK[i-1] = _IWORK(i);
    }
    for (INTEGER i=1; i<*N; ++i) {
        BWORK[i-1] = _BWORK(i);
    }
}

//-- zgeesx --------------------------------------------------------------------
void
LAPACK_DECL(zgeesx)(const char       *JOBVS,
                    const char       *SORT,
                    LOGICAL          (*SELECT)(const DOUBLE_COMPLEX *),
                    const char       *SENSE,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
                    INTEGER          *SDIM,
                    DOUBLE_COMPLEX   *W,
                    DOUBLE_COMPLEX   *VS,
                    const INTEGER    *LDVS,
                    DOUBLE           *RCONDE,
                    DOUBLE           *RCONDV,
                    DOUBLE_COMPLEX   *WORK,
                    const INTEGER    *LWORK,
                    DOUBLE           *RWORK,
                    LOGICAL          *BWORK,
                    INTEGER          *INFO)
{
    LAPACK_DEBUG_OUT("LAPACK INTERFACE: dgeesx");

    using std::max;
    using std::min;
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    const bool lQuery = (*LWORK==-1);
    const bool wantVS = (*JOBVS=='V');
    const bool wantST = (*SORT=='S');

    const bool wantSN = (*SENSE=='N');
    const bool wantSE = (*SENSE=='E');
    const bool wantSV = (*SENSE=='V');
    const bool wantSB = (*SENSE=='B');

    if ((!wantVS) && (*JOBVS!='N')) {
        *INFO = 1;
    } else if ((!wantST) && (*SORT!='N')) {
        *INFO = 2;
    } else if (!(wantSN || wantSE || wantSV || wantSB)
            || ( !wantST && !wantSN))
    {
        *INFO = 4;
    } else if (*N<0) {
        *INFO = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        *INFO = 7;
    } else if (*LDVS<1 || (wantVS && *LDVS<*N)) {
        *INFO = 11;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEESX", INFO);
        *INFO = -(*INFO);
        return;
    }

//
//  Setup FLENS matrix/vector types
//
    typedef DGeMatrixView::IndexType IndexType;

    SENSE::Sense sense = SENSE::Sense(*SENSE);

    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);
    auto zW     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(W);
    auto zVS    = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(VS);
    auto zWORK  = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    ZGeMatrixView       _A      = ZFSView(*N, *N, *LDA, zA);
    IndexType           _SDIM   = *SDIM;
    ZDenseVectorView    _W      = ZArrayView(*N, zW, 1);
    ZGeMatrixView       _VS     = ZFSView(*N, *N, *LDVS, zVS);
    ZDenseVectorView    _WORK   = ZArrayView(*LWORK, zWORK, 1);
    DDenseVectorView    _RWORK  = DArrayView(*N, RWORK, 1);

    BDenseVector        _BWORK(*N);
    for (INTEGER i=1; i<*N; ++i) {
        _BWORK(i) = BWORK[i-1];
    }

//
//  Test if work has at least minimal worksize
//
    auto wsQuery = esx_wsq(wantVS, sense, _A);
    INTEGER minWork = wsQuery.first;

    if (*LWORK<minWork && !lQuery) {
        *INFO = 15;
    }

    if (*INFO!=0) {
        LAPACK_ERROR("ZGEESX", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    SelectFunction<CXX_DOUBLE_COMPLEX> select(SELECT);

    esx(wantVS, wantST, select, sense, _A, _SDIM, _W, _VS,
        *RCONDE, *RCONDV, _WORK, _RWORK, _BWORK);

    *SDIM   = _SDIM;
    for (INTEGER i=1; i<*N; ++i) {
        BWORK[i-1] = _BWORK(i);
    }
}

} // extern "C"

} } // namespace lapack, flens
