#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *BETA,
            float           *Y,
            const INTEGER   *INCY)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _TRANS = toupper(*TRANS);

    if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 1;
    } else if (*M<0) {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*KL<0) {
        info = 4;
    } else if (*KU<0) {
        info = 5;
    } else if (*LDA<(*KL + *KU + 1)) {
        info = 8;
    } else if (*INCX==0) {
        info = 10;
    } else if (*INCY==0) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("SGBMV ", &info);
        return;
    }

    Transpose    trans = convertTo<Transpose>(_TRANS);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::gbmv(ColMajor, trans,
                  *M, *N,
                  *KL, *KU,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

void
BLAS(dgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const double    *ALPHA,
            const double    *_A,
            const INTEGER   *LDA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *BETA,
            double          *Y,
            const INTEGER   *INCY)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _TRANS = toupper(*TRANS);

    if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 1;
    } else if (*M<0) {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*KL<0) {
        info = 4;
    } else if (*KU<0) {
        info = 5;
    } else if (*LDA<(*KL + *KU + 1)) {
        info = 8;
    } else if (*INCX==0) {
        info = 10;
    } else if (*INCY==0) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("DGBMV ", &info);
        return;
    }

    Transpose    trans = convertTo<Transpose>(_TRANS);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::gbmv(ColMajor, trans,
                  *M, *N,
                  *KL, *KU,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

void
BLAS(cgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const cfloat    *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *BETA,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _TRANS = toupper(*TRANS);

    if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 1;
    } else if (*M<0) {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*KL<0) {
        info = 4;
    } else if (*KU<0) {
        info = 5;
    } else if (*LDA<(*KL + *KU + 1)) {
        info = 8;
    } else if (*INCX==0) {
        info = 10;
    } else if (*INCY==0) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("CGBMV ", &info);
        return;
    }

    Transpose    trans = convertTo<Transpose>(_TRANS);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::gbmv(ColMajor, trans,
                  *M, *N,
                  *KL, *KU,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

void
BLAS(zgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const cdouble   *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *BETA,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _TRANS = toupper(*TRANS);

    if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 1;
    } else if (*M<0) {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*KL<0) {
        info = 4;
    } else if (*KU<0) {
        info = 5;
    } else if (*LDA<(*KL + *KU + 1)) {
        info = 8;
    } else if (*INCX==0) {
        info = 10;
    } else if (*INCY==0) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("ZGBMV ", &info);
        return;
    }

    Transpose    trans = convertTo<Transpose>(_TRANS);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::gbmv(ColMajor, trans,
                  *M, *N,
                  *KL, *KU,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

} // extern "C"
