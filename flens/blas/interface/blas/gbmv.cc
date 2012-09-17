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
    const bool   noTrans = (trans==NoTrans || trans==Conj);
    INTEGER      lenX, lenY;

    if (noTrans) {
        lenX = *N;
        lenY = *M;
    } else {
        lenX = *M;
        lenY = *N;
    }

    SGbMatrixConstView    A = SBandConstView(*M, *N, *KL, *KU,  _A, *LDA);
    SDenseVectorConstView x(SConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
    SDenseVectorView      y(SArrayView(lenY, Y, abs(*INCY)), *INCY<0);

    blas::mv(trans, *ALPHA, A, x, *BETA, y);
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
    const bool   noTrans = (trans==NoTrans || trans==Conj);
    INTEGER      lenX, lenY;

    if (noTrans) {
        lenX = *N;
        lenY = *M;
    } else {
        lenX = *M;
        lenY = *N;
    }

    DGbMatrixConstView    A = DBandConstView(*M, *N, *KL, *KU,  _A, *LDA);
    DDenseVectorConstView x(DConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
    DDenseVectorView      y(DArrayView(lenY, Y, abs(*INCY)), *INCY<0);

    blas::mv(trans, *ALPHA, A, x, *BETA, y);
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
    const bool   noTrans = (trans==NoTrans || trans==Conj);
    INTEGER      lenX, lenY;

    if (noTrans) {
        lenX = *N;
        lenY = *M;
    } else {
        lenX = *M;
        lenY = *N;
    }

    CGbMatrixConstView    A = CBandConstView(*M, *N, *KL, *KU,  _A, *LDA);
    CDenseVectorConstView x(CConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
    CDenseVectorView      y(CArrayView(lenY, Y, abs(*INCY)), *INCY<0);

    blas::mv(trans, *ALPHA, A, x, *BETA, y);
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
    const bool   noTrans = (trans==NoTrans || trans==Conj);
    INTEGER      lenX, lenY;

    if (noTrans) {
        lenX = *N;
        lenY = *M;
    } else {
        lenX = *M;
        lenY = *N;
    }

    ZGbMatrixConstView    A = ZBandConstView(*M, *N, *KL, *KU,  _A, *LDA);
    ZDenseVectorConstView x(ZConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView      y(ZArrayView(lenY, Y, abs(*INCY)), *INCY<0);

    blas::mv(trans, *ALPHA, A, x, *BETA, y);
}

} // extern "C"
