#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sger)(const INTEGER   *M,
           const INTEGER   *N,
           const float     *ALPHA,
           const float     *X,
           const INTEGER   *INCX,
           const float     *Y,
           const INTEGER   *INCY,
           float           *A_,
           const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

        cblas_sger(CBLAS_ORDER::CblasColMajor,
                   *M, *N,
                   *ALPHA,
                   X, *INCX,
                   Y, *INCY,
                   A_, *LDA);
#   else

        using std::abs;
        using std::max;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;

            if (*M<0) {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*M)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("SGER  ", &info);
                return;
            }
#       endif

        SDenseVectorConstView x(SConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        SGeMatrixView         A = SFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(y);
#       else
            blas::r(*ALPHA, x, y, A);
#       endif
#   endif
}

void
BLAS(dger)(const INTEGER   *M,
           const INTEGER   *N,
           const double    *ALPHA,
           const double    *X,
           const INTEGER   *INCX,
           const double    *Y,
           const INTEGER   *INCY,
           double          *A_,
           const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

        cblas_dger(CBLAS_ORDER::CblasColMajor,
                   *M, *N,
                   *ALPHA,
                   X, *INCX,
                   Y, *INCY,
                   A_, *LDA);
#   else

        using std::abs;
        using std::max;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;

            if (*M<0) {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*M)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("DGER  ", &info);
                return;
            }
#       endif

        DDenseVectorConstView x(DConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        DDenseVectorConstView y(DConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        DGeMatrixView         A = DFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(y);
#       else
            blas::r(*ALPHA, x, y, A);
#       endif
#   endif

}

void
BLAS(cgerc)(const INTEGER   *M,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *Y,
            const INTEGER   *INCY,
            cfloat          *A_,
            const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_cgerc(CBLAS_ORDER::CblasColMajor,
                *M, *N,
                reinterpret_cast<const float *>(ALPHA),
                reinterpret_cast<const float *>(X), *INCX,
                reinterpret_cast<const float *>(Y), *INCY,
                reinterpret_cast<float *>(A_), *LDA);
#   else

        using std::abs;
        using std::max;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;

            if (*M<0) {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*M)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("CGERC ", &info);
                return;
            }
#       endif

        CDenseVectorConstView x(CConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        CGeMatrixView         A = CFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(y);
#       else
            blas::rc(*ALPHA, x, y, A);
#       endif
#   endif
}

void
BLAS(cgeru)(const INTEGER   *M,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *Y,
            const INTEGER   *INCY,
            cfloat          *A_,
            const INTEGER   *LDA)
{

#   ifdef TEST_DIRECT_CBLAS

        cblas_cgeru(CBLAS_ORDER::CblasColMajor,
                    *M, *N,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(X), *INCX,
                    reinterpret_cast<const float *>(Y), *INCY,
                    reinterpret_cast<float *>(A_), *LDA);
#   else

        using std::abs;
        using std::max;


#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;

            if (*M<0) {
                info = 1;
            } else if (*N<0) {
                info = 2;
            } else if (*INCX==0) {
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*M)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("CGERU ", &info);
                return;
            }
#       endif

        CDenseVectorConstView x(CConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        CDenseVectorConstView y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        CGeMatrixView         A = CFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(y);
#       else
            blas::ru(*ALPHA, x, y, A);
#       endif
#   endif

}

void
BLAS(zgerc)(const INTEGER   *M,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *Y,
            const INTEGER   *INCY,
            cdouble         *A_,
            const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_zgerc(CBLAS_ORDER::CblasColMajor,
                *M, *N,
                reinterpret_cast<const double *>(ALPHA),
                reinterpret_cast<const double *>(X), *INCX,
                reinterpret_cast<const double *>(Y), *INCY,
                reinterpret_cast<double *>(A_), *LDA);
#   else

        using std::abs;
        using std::max;

        INTEGER info  = 0;

        if (*M<0) {
            info = 1;
        } else if (*N<0) {
            info = 2;
        } else if (*INCX==0) {
            info = 5;
        } else if (*INCY==0) {
            info = 7;
        } else if (*LDA<max(INTEGER(1),*M)) {
            info = 9;
        }
        if (info!=0) {
            BLAS(xerbla)("ZGERC ", &info);
            return;
        }

        ZDenseVectorConstView x(ZConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        ZDenseVectorConstView y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        ZGeMatrixView         A = ZFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(y);
#       else
            blas::rc(*ALPHA, x, y, A);
#       endif
#   endif
}

void
BLAS(zgeru)(const INTEGER   *M,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *Y,
            const INTEGER   *INCY,
            cdouble         *A_,
            const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

    cblas_zgeru(CBLAS_ORDER::CblasColMajor,
                *M, *N,
                reinterpret_cast<const double *>(ALPHA),
                reinterpret_cast<const double *>(X), *INCX,
                reinterpret_cast<const double *>(Y), *INCY,
                reinterpret_cast<double *>(A_), *LDA);
#   else

        using std::abs;
        using std::max;

        INTEGER info  = 0;

        if (*M<0) {
            info = 1;
        } else if (*N<0) {
            info = 2;
        } else if (*INCX==0) {
            info = 5;
        } else if (*INCY==0) {
            info = 7;
        } else if (*LDA<max(INTEGER(1),*M)) {
            info = 9;
        }
        if (info!=0) {
            BLAS(xerbla)("ZGERU ", &info);
            return;
        }

        ZDenseVectorConstView x(ZConstArrayView(*M, X, abs(*INCX)), *INCX<0);
        ZDenseVectorConstView y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        ZGeMatrixView         A = ZFullView(*M, *N, A_, *LDA);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(y);
#       else
            blas::ru(*ALPHA, x, y, A);
#       endif
#   endif
}

} // extern "C"
