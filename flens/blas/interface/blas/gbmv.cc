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
            const float     *A_,
            const INTEGER   *LDA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *BETA,
            float           *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: sgbmv");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANS_ = toupper(*TRANS);
        const Transpose    trans  = convertTo<Transpose>(TRANS_);

        cblas_sgbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(trans),
                    *M, *N, *KL, *KU, *ALPHA, A_, *LDA, X, *INCX,
                    *BETA, Y, *INCY);

#   else

        using std::abs;
        using std::max;


        char    TRANS_ = toupper(*TRANS);

#       ifndef NO_INPUT_CHECK
            INTEGER info   = 0;
            if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
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
#       endif

        Transpose    trans = convertTo<Transpose>(TRANS_);
        const bool   noTrans = (trans==NoTrans || trans==Conj);
        INTEGER      lenX, lenY;

        if (noTrans) {
            lenX = *N;
            lenY = *M;
        } else {
            lenX = *M;
            lenY = *N;
        }

        SGbMatrixConstView    A = SBandConstView(*M, *N, *KL, *KU,  A_, *LDA);
        SDenseVectorConstView x(SConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
        SDenseVectorView      y(SArrayView(lenY, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (trans==NoTrans) {
                y = beta*y + alpha*A*x;
            } else if (trans==Trans) {
                y = beta*y + alpha*transpose(A)*x;
            } else if (trans==ConjTrans) {
                y = beta*y + alpha*conjTrans(A)*x;
            }
#       else
            blas::mv(trans, *ALPHA, A, x, *BETA, y);
#       endif
#   endif
}

void
BLAS(dgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *BETA,
            double          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dgbmv");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANS_ = toupper(*TRANS);
        const Transpose    trans  = convertTo<Transpose>(TRANS_);

        cblas_dgbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(trans),
                    *M, *N, *KL, *KU, *ALPHA, A_, *LDA, X, *INCX,
                    *BETA, Y, *INCY);


#   else

        using std::abs;
        using std::max;

        char    TRANS_ = toupper(*TRANS);

#       ifndef NO_INPUT_CHECK
            INTEGER info   = 0;
            if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
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
#       endif

        Transpose    trans = convertTo<Transpose>(TRANS_);
        const bool   noTrans = (trans==NoTrans || trans==Conj);
        INTEGER      lenX, lenY;

        if (noTrans) {
            lenX = *N;
            lenY = *M;
        } else {
            lenX = *M;
            lenY = *N;
        }

        DGbMatrixConstView    A = DBandConstView(*M, *N, *KL, *KU,  A_, *LDA);
        DDenseVectorConstView x(DConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
        DDenseVectorView      y(DArrayView(lenY, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (trans==NoTrans) {
                y = beta*y + alpha*A*x;
            } else if (trans==Trans) {
                y = beta*y + alpha*transpose(A)*x;
            } else if (trans==ConjTrans) {
                y = beta*y + alpha*conjTrans(A)*x;
            }
#       else
            blas::mv(trans, *ALPHA, A, x, *BETA, y);
#       endif
#   endif
}

void
BLAS(cgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *BETA,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cgbmv");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANS_ = toupper(*TRANS);
        const Transpose    trans  = convertTo<Transpose>(TRANS_);

        cblas_cgbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(trans),
                    *M, *N, *KL, *KU,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(X), *INCX,
                    reinterpret_cast<const float *>(BETA),
                    reinterpret_cast<float *>(Y), *INCY);

#   else

        using std::abs;
        using std::max;

        char    TRANS_ = toupper(*TRANS);

#       ifndef NO_INPUT_CHECK
            INTEGER info   = 0;
            if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
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
#       endif

        Transpose    trans = convertTo<Transpose>(TRANS_);
        const bool   noTrans = (trans==NoTrans || trans==Conj);
        INTEGER      lenX, lenY;

        if (noTrans) {
            lenX = *N;
            lenY = *M;
        } else {
            lenX = *M;
            lenY = *N;
        }

        CGbMatrixConstView    A = CBandConstView(*M, *N, *KL, *KU,  A_, *LDA);
        CDenseVectorConstView x(CConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
        CDenseVectorView      y(CArrayView(lenY, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (trans==NoTrans) {
                y = beta*y + alpha*A*x;
            } else if (trans==Trans) {
                y = beta*y + alpha*transpose(A)*x;
            } else if (trans==ConjTrans) {
                y = beta*y + alpha*conjTrans(A)*x;
            }
#       else
            blas::mv(trans, *ALPHA, A, x, *BETA, y);
#       endif
#   endif
}

void
BLAS(zgbmv)(const char      *TRANS,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *KL,
            const INTEGER   *KU,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *BETA,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zgbmv");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANS_ = toupper(*TRANS);
        const Transpose    trans  = convertTo<Transpose>(TRANS_);

        cblas_zgbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(trans),
                    *M, *N, *KL, *KU,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<const double *>(X), *INCX,
                    reinterpret_cast<const double *>(BETA),
                    reinterpret_cast<double *>(Y), *INCY);

#   else

        using std::abs;
        using std::max;

        char    TRANS_ = toupper(*TRANS);

#       ifndef NO_INPUT_CHECK
            INTEGER info   = 0;
            if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
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
#       endif
        Transpose    trans = convertTo<Transpose>(TRANS_);
        const bool   noTrans = (trans==NoTrans || trans==Conj);
        INTEGER      lenX, lenY;

        if (noTrans) {
            lenX = *N;
            lenY = *M;
        } else {
            lenX = *M;
            lenY = *N;
        }

        ZGbMatrixConstView    A = ZBandConstView(*M, *N, *KL, *KU,  A_, *LDA);
        ZDenseVectorConstView x(ZConstArrayView(lenX, X, abs(*INCX)), *INCX<0);
        ZDenseVectorView      y(ZArrayView(lenY, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (trans==NoTrans) {
                y = beta*y + alpha*A*x;
            } else if (trans==Trans) {
                y = beta*y + alpha*transpose(A)*x;
            } else if (trans==ConjTrans) {
                y = beta*y + alpha*conjTrans(A)*x;
            }
#       else
            blas::mv(trans, *ALPHA, A, x, *BETA, y);
#       endif
#   endif
}

} // extern "C"
