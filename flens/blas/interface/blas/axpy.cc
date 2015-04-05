#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(saxpy)(const INTEGER   *N,
            const float     *ALPHA,
            const float     *X,
            const INTEGER   *INCX,
            float           *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: saxpy");

    using std::abs;

#   ifdef TEST_DIRECT_CBLAS

        cblas_saxpy(*N, *ALPHA, X, *INCX, Y, *INCY);

#   else

        SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorView       y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            y += (*ALPHA) * x;
#       else
            blas::axpy(*ALPHA, x, y);
#       endif
#   endif
}

void
BLAS(daxpy)(const INTEGER   *N,
            const double    *ALPHA,
            const double    *X,
            const INTEGER   *INCX,
            double          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: daxpy");

    using std::abs;

#   ifdef TEST_DIRECT_CBLAS

    cblas_daxpy(*N, *ALPHA, X, *INCX, Y, *INCY);

#   else

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorView       y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            y += (*ALPHA) * x;
#       else
            blas::axpy(*ALPHA, x, y);
#       endif
#   endif
}

void
BLAS(caxpy)(const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *X,
            const INTEGER   *INCX,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: caxpy");

    using std::abs;

#   ifdef TEST_DIRECT_CBLAS

    cblas_caxpy(*N, reinterpret_cast<const float *>(ALPHA),
                reinterpret_cast<const float *>(X), *INCX,
                reinterpret_cast<float *>(Y), *INCY);

#   else

        CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CDenseVectorView       y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            y += (*ALPHA) * x;
#       else
            blas::axpy(*ALPHA, x, y);
#       endif
#   endif
}

void
BLAS(zaxpy)(const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *X,
            const INTEGER   *INCX,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zaxpy");

    using std::abs;

#   ifdef TEST_DIRECT_CBLAS

    cblas_zaxpy(*N, reinterpret_cast<const double *>(ALPHA),
                reinterpret_cast<const double *>(X), *INCX,
                reinterpret_cast<double *>(Y), *INCY);

#   else

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZDenseVectorView       y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            y += (*ALPHA) * x;
#       else
            blas::axpy(*ALPHA, x, y);
#       endif
#   endif
}

} // extern "C"
